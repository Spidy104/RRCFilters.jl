"""
    carriersync(x; modulation=:qam, sps=2, damping_factor=0.707, loop_bandwidth=0.01)

Recover carrier phase/frequency offset with a hard-decision digital PLL
(Costas loop). Returns `(corrected, phase_estimate)`:
`corrected` is always complex regardless of whether `x` is real or complex,
and `phase_estimate` is the running phase estimate in radians.

`modulation` selects the phase-error detector: `:qam` (any square/cross QAM
order, or any 4-point constellation aligned to the diagonals, e.g.
`qammod(bits, 4)` or `pskmod(bits, 4; phase_offset=pi/4)`), `:bpsk`
(2-point, axis-aligned, e.g. `pskmod(bits, 2)`), or `:psk8` (8-point,
axis-and-diagonal-aligned, e.g. `pskmod(bits, 8)`). `sps` is
`SamplesPerSymbol`, a loop-gain scaling parameter -- `carriersync` does not
resample; pass `sps=1` for an already symbol-rate stream. The internal loop
state always starts at zero.
An empty `x` returns two empty vectors.
"""
function carriersync end

function carriersync(x::AbstractVector{<:Number}; modulation::Symbol=:qam, sps::Integer=2,
                      damping_factor::Real=0.707, loop_bandwidth::Real=0.01)
    _validate_carriersync_arguments(modulation, sps, damping_factor, loop_bandwidth)

    if isempty(x)
        real_type = _awgn_real_type(_upfirdn_empty_type(eltype(x)))
        return Complex{real_type}[], real_type[]
    end

    output_type = _awgn_output_type(x)
    real_type = _awgn_real_type(output_type)

    kp = _carriersync_modulation_gain(Val(modulation), real_type)
    gains = _carriersync_loop_gains(real_type(damping_factor), real_type(loop_bandwidth), sps, kp)

    corrected = Vector{Complex{real_type}}(undef, length(x))
    phase_estimate = Vector{real_type}(undef, length(x))
    _carriersync_kernel!(corrected, phase_estimate, x, gains, Val(modulation))
    return corrected, phase_estimate
end

function _validate_carriersync_arguments(modulation::Symbol, sps::Integer, damping_factor::Real, loop_bandwidth::Real)
    modulation in (:qam, :bpsk, :psk8) || throw(ArgumentError("modulation must be :qam, :bpsk, or :psk8"))

    sps isa Bool && throw(ArgumentError("sps must be an integer, not Bool"))
    sps > 0 || throw(ArgumentError("sps must be positive"))

    damping_factor isa Bool && throw(ArgumentError("damping_factor must be a real number, not Bool"))
    isfinite(damping_factor) || throw(ArgumentError("damping_factor must be finite"))
    damping_factor > 0 || throw(ArgumentError("damping_factor must be positive"))

    loop_bandwidth isa Bool && throw(ArgumentError("loop_bandwidth must be a real number, not Bool"))
    isfinite(loop_bandwidth) || throw(ArgumentError("loop_bandwidth must be finite"))
    (loop_bandwidth > 0 && loop_bandwidth <= 1) || throw(ArgumentError("loop_bandwidth must be in (0, 1]"))

    return nothing
end

_carriersync_modulation_gain(::Val{:qam}, ::Type{T}) where {T} = T(2)
_carriersync_modulation_gain(::Val{:bpsk}, ::Type{T}) where {T} = T(1)
_carriersync_modulation_gain(::Val{:psk8}, ::Type{T}) where {T} = T(1)

const _CARRIERSYNC_INV_SQRT2 = inv(sqrt(2.0))
const _CARRIERSYNC_PSK8_SECTOR_RATIO = sqrt(2.0) - 1

# Rice's second-order proportional-integrator digital PLL. Theta is computed
# via this literal two-step form (Bn*K0, then divided by (zeta+0.25/zeta)*K0)
# than the algebraically-simplified Bn/(zeta+0.25/zeta), and both gains
# DIVIDE by Kp*K0 -- not multiply, which would give qualitatively wrong loop
# dynamics, not just a rounding difference.
function _carriersync_loop_gains(damping_factor::T, loop_bandwidth::T, sps::Integer, kp::T) where {T<:AbstractFloat}
    k0 = T(sps)
    theta = (loop_bandwidth * k0) / ((damping_factor + T(0.25) / damping_factor) * k0)
    d = 1 + 2 * damping_factor * theta + theta * theta
    proportional_gain = (4 * damping_factor * theta / d) / (kp * k0)
    integrator_gain = (4 * theta * theta / d) / (kp * k0)
    return (proportional_gain=proportional_gain, integrator_gain=integrator_gain)
end

_carriersync_phase_error(s::Complex{T}, ::Val{:qam}) where {T} = sign(real(s)) * imag(s) - sign(imag(s)) * real(s)
_carriersync_phase_error(s::Complex{T}, ::Val{:bpsk}) where {T} = sign(real(s)) * imag(s)

function _carriersync_phase_error(s::Complex{T}, ::Val{:psk8}) where {T}
    re, im = real(s), imag(s)
    are, aim = abs(re), abs(im)
    ratio = T(_CARRIERSYNC_PSK8_SECTOR_RATIO)
    invsqrt2 = T(_CARRIERSYNC_INV_SQRT2)
    if are >= aim
        if aim < ratio * are
            return re >= zero(T) ? im : -im
        end
        return re >= zero(T) ? (im >= zero(T) ? invsqrt2 * (im - re) : invsqrt2 * (im + re)) :
                               (im >= zero(T) ? -invsqrt2 * (im + re) : invsqrt2 * (re - im))
    end
    if are < ratio * aim
        return im >= zero(T) ? -re : re
    end
    return re >= zero(T) ? (im >= zero(T) ? invsqrt2 * (im - re) : invsqrt2 * (im + re)) :
                           (im >= zero(T) ? -invsqrt2 * (im + re) : invsqrt2 * (re - im))
end

# In the per-sample recursion, both `phase` and `previous_sample` used inside
# iteration k are the values
# produced at the END of iteration k-1 (or zero initial state for k=1) --
# `dds_acc` lags `dds_input` by one further iteration, so a phase-error
# computed this sample doesn't reach the DDS phase accumulator until the next
# one. `output`/`phase_estimate` are freshly allocated 1-based vectors, so
# `enumerate(x)` (which always yields a 1-based index regardless of `x`'s own
# indexing) is the correct generic-over-indexing pattern here.
function _carriersync_kernel!(output::Vector{Complex{T}}, phase_estimate::Vector{T},
                               x::AbstractVector{<:Number}, gains, modtag::Val) where {T<:AbstractFloat}
    phase = zero(T)
    previous_sample = zero(Complex{T})
    integral_acc = zero(T)
    dds_acc = zero(T)
    dds_input = zero(T)
    for (index, sample) in enumerate(x)
        pherr = _carriersync_phase_error(previous_sample, modtag)
        corrected = Complex{T}(sample) * cis(phase)
        output[index] = corrected
        integral_acc = pherr * gains.integrator_gain + integral_acc
        dds_acc = dds_input + dds_acc
        dds_input = pherr * gains.proportional_gain + integral_acc
        phase = -dds_acc
        phase_estimate[index] = dds_acc
        previous_sample = corrected
    end
    return output, phase_estimate
end
