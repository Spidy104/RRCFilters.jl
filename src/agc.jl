"""
    agc(x; desired_amplitude=1.0, step_size=0.01, initial_gain=1.0, min_gain=1e-4, max_gain=1e4)

Normalize a signal's amplitude toward `desired_amplitude` with a feedback-
driven Automatic Gain Control loop, the classical receiver front-end stage
that runs before carrier/timing synchronization. Returns `(y, gain)`:
`y[i] = gain[i] * x[i]` for every sample (no stuffing/skipping, unlike
[`symbolsync`](@ref) -- output length is always exactly `length(x)`), and
`gain[i]` is the nonnegative real scalar gain actually applied to `x[i]`,
computed from state available strictly before that sample. Realness is
preserved exactly as in [`timingoffset`](@ref): real `x` in gives real `y`
out, complex in gives complex out.

This is a transparent first-principles design with explicit loop state and
no hidden persistent object behavior.

The loop tracks the error between `desired_amplitude` and the current
sample's magnitude in the LOG domain, with `step_size` (required in
`(0, 1]` for a non-oscillating single-pole response; `step_size=1` tracks
the instantaneous target with no smoothing at all) as the pole's rate.
Operating in the log domain, rather than a naive linear-domain error loop,
makes the loop's convergence time (in samples) independent of the input's
absolute amplitude -- a real, textbook AGC property. `gain` is clamped to
`[min_gain, max_gain]` every sample, in the log domain, before
exponentiating, so the tracked state itself can never overflow;
`initial_gain` (which must lie in `[min_gain, max_gain]`) is the gain
applied to the very first sample, before any adaptation has happened. A
sample whose magnitude is extreme enough that even a clamped gain would
overflow the output's floating-point range is handled by transparently
reducing the gain actually applied to that one sample (and reporting that
reduced value in `gain`, so `y[i] == gain[i] * x[i]` always holds exactly)
-- the loop's internal state is not corrupted by this, since the state
update always uses the true gain trajectory, not the realized output.

An all-zero (silent) input drives `gain` to saturate at `max_gain` within
a few samples, never to `Inf` or `NaN`. An empty `x` returns two empty
vectors.
"""
function agc end

function agc(x::AbstractVector{<:Number}; desired_amplitude::Real=1.0, step_size::Real=0.01,
             initial_gain::Real=1.0, min_gain::Real=1e-4, max_gain::Real=1e4)
    _validate_agc_arguments(desired_amplitude, step_size, initial_gain, min_gain, max_gain)

    if isempty(x)
        output_type = _upfirdn_empty_type(eltype(x))
        real_type = _awgn_real_type(output_type)
        return output_type[], real_type[]
    end

    output_type = _awgn_output_type(x)
    real_type = _awgn_real_type(output_type)
    constants = _agc_loop_constants(real_type(desired_amplitude), real_type(step_size),
                                     real_type(initial_gain), real_type(min_gain), real_type(max_gain))

    y = Vector{output_type}(undef, length(x))
    gain = Vector{real_type}(undef, length(x))
    _agc_kernel!(y, gain, x, constants)
    return y, gain
end

function _validate_agc_arguments(desired_amplitude::Real, step_size::Real, initial_gain::Real,
                                  min_gain::Real, max_gain::Real)
    desired_amplitude isa Bool && throw(ArgumentError("desired_amplitude must be a real number, not Bool"))
    isfinite(desired_amplitude) || throw(ArgumentError("desired_amplitude must be finite"))
    desired_amplitude > 0 || throw(ArgumentError("desired_amplitude must be positive"))

    step_size isa Bool && throw(ArgumentError("step_size must be a real number, not Bool"))
    isfinite(step_size) || throw(ArgumentError("step_size must be finite"))
    (step_size > 0 && step_size <= 1) || throw(ArgumentError("step_size must be in (0, 1]"))

    initial_gain isa Bool && throw(ArgumentError("initial_gain must be a real number, not Bool"))
    isfinite(initial_gain) || throw(ArgumentError("initial_gain must be finite"))
    initial_gain > 0 || throw(ArgumentError("initial_gain must be positive"))

    min_gain isa Bool && throw(ArgumentError("min_gain must be a real number, not Bool"))
    isfinite(min_gain) || throw(ArgumentError("min_gain must be finite"))
    min_gain > 0 || throw(ArgumentError("min_gain must be positive"))

    max_gain isa Bool && throw(ArgumentError("max_gain must be a real number, not Bool"))
    isfinite(max_gain) || throw(ArgumentError("max_gain must be finite"))
    max_gain > 0 || throw(ArgumentError("max_gain must be positive"))

    min_gain <= max_gain || throw(ArgumentError("min_gain must be <= max_gain"))
    (min_gain <= initial_gain <= max_gain) || throw(ArgumentError("initial_gain must be within [min_gain, max_gain]"))
    return nothing
end

function _agc_loop_constants(desired_amplitude::T, step_size::T, initial_gain::T,
                              min_gain::T, max_gain::T) where {T<:AbstractFloat}
    log_min_gain = log(min_gain)
    log_max_gain = log(max_gain)
    log_g0 = clamp(log(initial_gain), log_min_gain, log_max_gain)
    log_desired = log(desired_amplitude)
    log_fmax_margin = log(floatmax(T)) - log(T(4))
    return (log_desired=log_desired, step_size=step_size, log_min_gain=log_min_gain,
            log_max_gain=log_max_gain, log_g0=log_g0, log_fmax_margin=log_fmax_margin)
end

function _agc_kernel!(y::Vector{Ty}, gain::Vector{T}, x::AbstractVector{<:Number}, c) where {Ty<:Number,T<:AbstractFloat}
    log_g = c.log_g0
    fmax = floatmax(T)
    fmin = floatmin(T)
    for (index, sample) in enumerate(x)
        converted = Ty(sample)
        mag_x = T(abs(converted))
        safe_mag_x = max(mag_x, fmin)
        log_mag_x = log(safe_mag_x)

        log_g_applied = min(log_g, c.log_fmax_margin - log_mag_x)
        g_applied = exp(log_g_applied)
        candidate = g_applied * converted
        if !isfinite(candidate)
            g_applied = fmax / safe_mag_x
            candidate = g_applied * converted
        end
        y[index] = candidate
        gain[index] = g_applied

        e = c.log_desired - log_g - log_mag_x
        log_g = clamp(log_g + c.step_size * e, c.log_min_gain, c.log_max_gain)
    end
    return y, gain
end
