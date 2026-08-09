"""
    symbolsync(x; sps, timing_error_detector=:zero_crossing, damping_factor=1.0,
               loop_bandwidth=0.01, detector_gain=2.7)

Recover symbol timing with a Farrow-interpolator PLL. Returns
`(symbols, timing_error)`: `symbols` is
always complex regardless of whether `x` is real or complex, and its
length is inherently variable -- nominally `length(x)/sps`, but the loop's
buffer stuffing/skipping mechanism can add or drop a symbol while tracking a
drifting sample-to-symbol clock ratio.
`timing_error` has length `length(x)` exactly (one value per input sample,
not per output symbol) and holds the running fractional timing estimate,
in `[0, 1)`.

`sps` is the required samples-per-symbol ratio.
`timing_error_detector` selects the error detector: `:zero_crossing`
(decision-directed, the default), `:gardner` (non-data-aided),
`:early_late` (non-data-aided), or `:mueller_muller` (decision-directed).
`damping_factor`, `loop_bandwidth`, and `detector_gain` configure the loop
filter. The internal loop state always starts at zero. An empty `x` returns two
empty vectors.

The inverse operation is [`timingoffset`](@ref), which introduces (rather
than recovers) a fractional sample-timing offset.
"""
function symbolsync end

function symbolsync(x::AbstractVector{<:Number}; sps::Integer,
                     timing_error_detector::Symbol=:zero_crossing,
                     damping_factor::Real=1.0, loop_bandwidth::Real=0.01,
                     detector_gain::Real=2.7)
    _validate_symbolsync_arguments(sps, timing_error_detector, damping_factor, loop_bandwidth, detector_gain)

    if isempty(x)
        real_type = _awgn_real_type(_upfirdn_empty_type(eltype(x)))
        return Complex{real_type}[], real_type[]
    end

    output_type = _awgn_output_type(x)
    real_type = _awgn_real_type(output_type)

    gains = _symbolsync_loop_gains(real_type(damping_factor), real_type(loop_bandwidth),
                                    real_type(detector_gain), sps)

    symbols = Vector{Complex{real_type}}(undef, 0)
    sizehint!(symbols, cld(length(x) * 11, 10 * sps))
    timing_error = Vector{real_type}(undef, length(x))
    _symbolsync_kernel!(symbols, timing_error, x, Int(sps), gains, Val(timing_error_detector))
    return symbols, timing_error
end

function _validate_symbolsync_arguments(sps::Integer, timing_error_detector::Symbol, damping_factor::Real,
                                         loop_bandwidth::Real, detector_gain::Real)
    timing_error_detector in (:zero_crossing, :gardner, :early_late, :mueller_muller) ||
        throw(ArgumentError("timing_error_detector must be :zero_crossing, :gardner, :early_late, or :mueller_muller"))

    sps isa Bool && throw(ArgumentError("sps must be an integer, not Bool"))
    sps >= 2 || throw(ArgumentError("sps must be >= 2"))

    damping_factor isa Bool && throw(ArgumentError("damping_factor must be a real number, not Bool"))
    isfinite(damping_factor) || throw(ArgumentError("damping_factor must be finite"))
    damping_factor > 0 || throw(ArgumentError("damping_factor must be positive"))

    loop_bandwidth isa Bool && throw(ArgumentError("loop_bandwidth must be a real number, not Bool"))
    isfinite(loop_bandwidth) || throw(ArgumentError("loop_bandwidth must be finite"))
    (loop_bandwidth > 0 && loop_bandwidth < 1) || throw(ArgumentError("loop_bandwidth must be in (0, 1)"))

    detector_gain isa Bool && throw(ArgumentError("detector_gain must be a real number, not Bool"))
    isfinite(detector_gain) || throw(ArgumentError("detector_gain must be finite"))
    detector_gain > 0 || throw(ArgumentError("detector_gain must be positive"))

    return nothing
end

function _symbolsync_loop_gains(damping_factor::T, loop_bandwidth::T, detector_gain::T, sps::Integer) where {T<:AbstractFloat}
    k0 = T(-1)
    theta = loop_bandwidth / T(sps) / (damping_factor + T(0.25) / damping_factor)
    d = (1 + 2 * damping_factor * theta + theta * theta) * k0 * detector_gain
    proportional_gain = (4 * damping_factor * theta) / d
    integrator_gain = (4 * theta * theta) / d
    return (proportional_gain=proportional_gain, integrator_gain=integrator_gain)
end

_symbolsync_buffer_indices(sps::Integer) = (fld(sps, 2) + 1 - rem(sps, 2), fld(sps, 2) + 1)


_symbolsync_ted_error(buf, mids, x, ::Val{:zero_crossing}) =
    real(mids) * (sign(real(buf[1])) - sign(real(x))) + imag(mids) * (sign(imag(buf[1])) - sign(imag(x)))
_symbolsync_ted_error(buf, mids, x, ::Val{:gardner}) =
    real(mids) * (real(buf[1]) - real(x)) + imag(mids) * (imag(buf[1]) - imag(x))
_symbolsync_ted_error(buf, mids, x, ::Val{:early_late}) =
    real(buf[end]) * (real(x) - real(buf[end-1])) + imag(buf[end]) * (imag(x) - imag(buf[end-1]))
_symbolsync_ted_error(buf, mids, x, ::Val{:mueller_muller}) =
    sign(real(buf[1])) * real(x) - sign(real(x)) * real(buf[1]) + sign(imag(buf[1])) * imag(x) - sign(imag(x)) * imag(buf[1])

_symbolsync_trigger(strobe, history, ::Val{:zero_crossing}) = strobe && !any(@view history[2:end])
_symbolsync_trigger(strobe, history, ::Val{:gardner}) = strobe && !any(@view history[2:end])
_symbolsync_trigger(strobe, history, ::Val{:mueller_muller}) = strobe && !any(@view history[2:end])
_symbolsync_trigger(strobe, history, ::Val{:early_late}) = !strobe && history[end] && !any(@view history[1:end-1])

function _symbolsync_shift1!(buffer::Vector{V}, value::V) where {V}
    n = length(buffer)
    copyto!(buffer, 1, buffer, 2, n - 1)
    buffer[n] = value
    return buffer
end

function _symbolsync_shift2!(buffer::Vector{Complex{T}}, value::Complex{T}) where {T}
    n = length(buffer)
    n > 2 && copyto!(buffer, 1, buffer, 3, n - 2)
    buffer[n-1] = zero(Complex{T})
    buffer[n] = value
    return buffer
end

function _symbolsync_kernel!(symbols::Vector{Complex{T}}, timing_error::Vector{T},
                              x::AbstractVector{<:Number}, sps::Int, gains, tedtag::Val) where {T<:AbstractFloat}
    mid, mid2 = _symbolsync_buffer_indices(sps)
    ted_buffer = zeros(Complex{T}, sps)
    history = fill(false, sps)
    s1 = s2 = s3 = zero(Complex{T})
    loop_state = zero(T)
    mu = zero(T)
    nco = zero(T)
    strobe = false
    inv_sps = one(T) / T(sps)

    for (index, sample) in enumerate(x)
        timing_error[index] = mu
        xc = Complex{T}(sample)
        y, s1, s2, s3 = _farrow_interp(xc, s1, s2, s3, mu)
        strobe && push!(symbols, y)

        mid_sample = (ted_buffer[mid] + ted_buffer[mid2]) / 2
        triggered = _symbolsync_trigger(strobe, history, tedtag)
        e = triggered ? _symbolsync_ted_error(ted_buffer, mid_sample, y, tedtag) : zero(T)

        strobe_count = Int(strobe) + count(identity, @view history[2:end])
        strobe_count == 1 && _symbolsync_shift1!(ted_buffer, y)
        strobe_count > 1 && _symbolsync_shift2!(ted_buffer, y)

        loop_out = e * gains.integrator_gain + loop_state
        v = e * gains.proportional_gain + loop_out
        loop_state = loop_out

        w = v + inv_sps
        _symbolsync_shift1!(history, strobe)
        new_strobe = nco < w
        new_strobe && (mu = nco / w)
        nco = mod(nco - w, one(T))
        strobe = new_strobe
    end
    return symbols, timing_error
end
