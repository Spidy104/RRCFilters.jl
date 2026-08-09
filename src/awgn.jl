"""
    awgn(x, snr_db; signal_power=nothing)
    awgn(rng, x, snr_db; signal_power=nothing)

Add white Gaussian noise to a real or complex signal at the requested SNR in
dB. By default, signal RMS is measured from the samples; `signal_power`
supplies a nonnegative linear-power override. Empty signals, non-finite values,
and values not representable in the selected output precision are rejected.
The result is a newly allocated vector and does not modify `x`. Concrete
floating-point inputs preserve their precision; integer and rational inputs
promote to `Float64`. For a zero measured or overridden power, this function
returns a noiseless copy because SNR is undefined. The overload accepting
`rng::AbstractRNG` permits reproducible simulations with RNGs that support
`randn(rng, T, dims)`. Real and imaginary complex-noise components have equal
power.
"""
function awgn end

function awgn(x::AbstractVector{<:Number}, snr_db::Real; signal_power::Union{Real,Nothing}=nothing)
    return awgn(default_rng(), x, snr_db; signal_power)
end

function awgn(rng::AbstractRNG, x::AbstractVector{T}, snr_db::Real; signal_power::Union{Real,Nothing}=nothing) where {T<:Number}
    _validate_awgn_inputs(x, snr_db, signal_power)

    output_type = _awgn_output_type(x)
    real_type = _awgn_real_type(output_type)
    work_type = _awgn_work_type(real_type, signal_power)
    signal_rms = isnothing(signal_power) ? _awgn_rms(x, output_type, work_type) : sqrt(work_type(signal_power))
    isfinite(signal_rms) || throw(ArgumentError("signal RMS must be finite"))

    snr = work_type(snr_db)
    isfinite(snr) || throw(ArgumentError("snr_db is not representable in the working precision"))
    scale_work = signal_rms * exp10(-snr / work_type(20))
    isfinite(scale_work) || throw(ArgumentError("requested noise amplitude is not finite"))
    scale = real_type(scale_work)
    isfinite(scale) || throw(ArgumentError("requested noise amplitude is not representable in the output precision"))
    iszero(scale) && return output_type.(x)

    noisy = randn(rng, output_type, length(x))
    for (index, sample) in enumerate(x)
        noisy[index] = output_type(sample) + scale * noisy[index]
    end
    return noisy
end

function _awgn_work_type(real_type::Type{<:AbstractFloat}, signal_power)
    return real_type === BigFloat || (!isnothing(signal_power) && !(signal_power isa Union{Float16,Float32,Float64})) ? BigFloat : Float64
end

function _awgn_rms(x::AbstractVector{<:Number}, output_type::Type, work_type::Type{<:AbstractFloat})
    fast = _awgn_rms_fast(x, output_type, work_type)
    isnothing(fast) || return fast
    return _awgn_rms_safe(x, output_type, work_type)
end

function _awgn_rms_fast(x::AbstractVector{<:Number}, output_type::Type, work_type::Type{<:AbstractFloat})
    n = length(x)
    x_first = firstindex(x)
    sum_of_squares = zero(work_type)
    valid = true
    any_nonzero = false
    @inbounds @simd for offset in 0:(n - 1)
        sample = x[x_first + offset]
        converted = output_type(sample)
        magnitude = work_type(abs(converted))
        valid &= isfinite(magnitude)
        valid &= !(!iszero(sample) && iszero(converted))
        any_nonzero |= !iszero(magnitude)
        sum_of_squares += magnitude * magnitude
    end
    (valid && isfinite(sum_of_squares)) || return nothing
    # A nonzero magnitude squaring to exactly zero (e.g. 1e-200) is the same
    # catastrophic-underflow failure mode _awgn_rms_safe's online rescale
    # exists to avoid, just at the opposite end of the range from overflow;
    # this path never squares a raw magnitude, so it can't happen there.
    (any_nonzero && iszero(sum_of_squares)) && return nothing
    return sqrt(sum_of_squares / work_type(n))
end

function _awgn_rms_safe(x::AbstractVector{<:Number}, output_type::Type, work_type::Type{<:AbstractFloat})
    scale = zero(work_type)
    sum_of_squares = one(work_type)

    for sample in x
        converted = output_type(sample)
        isfinite(converted) || throw(ArgumentError("x contains values not representable in the output precision"))
        !iszero(sample) && iszero(converted) && throw(ArgumentError("x contains values too small for the output precision"))
        magnitude = work_type(abs(converted))
        iszero(magnitude) && continue

        if scale < magnitude
            sum_of_squares = one(work_type) + sum_of_squares * (scale / magnitude)^2
            scale = magnitude
        else
            sum_of_squares += (magnitude / scale)^2
        end
    end

    return scale * sqrt(sum_of_squares / work_type(length(x)))
end

function _validate_awgn_inputs(x::AbstractVector{<:Number}, snr_db::Real, signal_power::Union{Real,Nothing})
    isempty(x) && throw(ArgumentError("x must be nonempty"))
    snr_db isa Bool && throw(ArgumentError("snr_db must be a real number, not Bool"))
    isfinite(snr_db) || throw(ArgumentError("snr_db must be finite"))

    if !isnothing(signal_power)
        signal_power isa Bool && throw(ArgumentError("signal_power must be a real number, not Bool"))
        isfinite(signal_power) || throw(ArgumentError("signal_power must be finite"))
        signal_power >= zero(signal_power) || throw(ArgumentError("signal_power must be nonnegative"))
    end
end
