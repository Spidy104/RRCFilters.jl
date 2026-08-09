"""
    freqoffset(x; freq_offset=0.0, phase_offset=0.0, sample_rate=1.0)

Apply a carrier frequency and phase offset impairment. The output is
always complex regardless of whether `x` is real or complex, since rotating
a real signal off the real axis is the function's whole purpose.

`freq_offset` and `sample_rate` are in Hz (`sample_rate` defaults to 1, so
`freq_offset` doubles as cycles-per-sample at the default). `phase_offset`
is in radians, keeping every angle in this library in one unit
(`pskmod`/`pskdemod`'s `phase_offset`, `carriersync`'s `phase_estimate`).
An empty `x` returns an empty complex vector.

The inverse operation is [`carriersync`](@ref), which recovers (rather than
applies) a phase/frequency offset.
"""
function freqoffset end

function freqoffset(x::AbstractVector{<:Number}; freq_offset::Real=0.0,
                     phase_offset::Real=0.0, sample_rate::Real=1.0)
    _validate_freqoffset_arguments(freq_offset, phase_offset, sample_rate)

    if isempty(x)
        real_type = _awgn_real_type(_upfirdn_empty_type(eltype(x)))
        return Complex{real_type}[]
    end

    output_type = _awgn_output_type(x)
    real_type = _awgn_real_type(output_type)

    frequency = real_type(freq_offset)
    rate = real_type(sample_rate)
    phase = real_type(phase_offset)
    isfinite(frequency) || throw(ArgumentError("freq_offset is not representable in the output precision"))
    isfinite(rate) || throw(ArgumentError("sample_rate is not representable in the output precision"))
    isfinite(phase) || throw(ArgumentError("phase_offset is not representable in the output precision"))

    rf, rfs = _freqoffset_rat(rem(frequency, rate) / rate, eps(real_type))
    phase_offset_pi = phase / real_type(pi)

    y = Vector{Complex{real_type}}(undef, length(x))
    _freqoffset_kernel!(y, x, rf, rfs, phase_offset_pi)
    return y
end

function _validate_freqoffset_arguments(freq_offset::Real, phase_offset::Real, sample_rate::Real)
    freq_offset isa Bool && throw(ArgumentError("freq_offset must be a real number, not Bool"))
    isfinite(freq_offset) || throw(ArgumentError("freq_offset must be finite"))

    phase_offset isa Bool && throw(ArgumentError("phase_offset must be a real number, not Bool"))
    isfinite(phase_offset) || throw(ArgumentError("phase_offset must be finite"))

    sample_rate isa Bool && throw(ArgumentError("sample_rate must be a real number, not Bool"))
    isfinite(sample_rate) || throw(ArgumentError("sample_rate must be finite"))
    sample_rate > 0 || throw(ArgumentError("sample_rate must be positive"))

    return nothing
end

function _freqoffset_rat(target::T, tol::T) where {T<:AbstractFloat}
    isfinite(target) || throw(ArgumentError("frequency ratio must be finite"))
    isfinite(tol) && tol > zero(T) || throw(ArgumentError("frequency-ratio tolerance must be positive and finite"))
    n1, n0 = one(T), zero(T)
    d1, d0 = zero(T), one(T)
    val = target
    while true
        d = round(val, RoundNearestTiesAway)
        val -= d
        n1, n0 = d * n1 + n0, n1
        d1, d0 = d * d1 + d0, d1
        (iszero(val) || abs(n1 / d1 - target) <= max(tol, eps(target))) && break
        val = 1 / val
    end
    return n1 / sign(d1), abs(d1)
end

const _FREQOFFSET_TABLE_LIMIT = 62500

function _freqoffset_kernel!(y::Vector{Complex{T}}, x::AbstractVector{<:Number},
                              rf::T, rfs::T, phase_offset_pi::T) where {T<:AbstractFloat}
    if rfs <= _FREQOFFSET_TABLE_LIMIT && rfs <= length(x)
        _freqoffset_table_kernel!(y, x, Int(mod(rf, rfs)), Int(rfs), rfs, phase_offset_pi)
    else
        for (index, sample) in enumerate(x)
            idx = mod(rf * T(index - 1), rfs)
            s, c = sincospi(2 * idx / rfs + phase_offset_pi)
            y[index] = Complex{T}(sample) * Complex{T}(c, s)
        end
    end
    return y
end

function _freqoffset_table_kernel!(y::Vector{Complex{T}}, x::AbstractVector{<:Number},
                                    step::Int, rfs_int::Int, rfs::T, phase_offset_pi::T) where {T<:AbstractFloat}
    table = Vector{Complex{T}}(undef, rfs_int)
    for k in 0:(rfs_int - 1)
        s, c = sincospi(2 * T(k) / rfs + phase_offset_pi)
        table[k + 1] = Complex{T}(c, s)
    end
    idx = 0
    @inbounds for (index, sample) in enumerate(x)
        y[index] = Complex{T}(sample) * table[idx + 1]
        idx += step
        idx >= rfs_int && (idx -= rfs_int)
    end
    return y
end
