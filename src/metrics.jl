"""
    biterr(tx_bits, rx_bits; mode=:bit, M=nothing)

Return `(error_count, error_rate)` for equal-length, nonempty bit vectors.
Set `mode=:symbol` and provide a power-of-two `M` to count a group with one or
more incorrect bits as one symbol error.
"""
function biterr end

"""
    evm(ideal_symbols, rx_symbols; mode=:percent)

Return normalized RMS error-vector magnitude as a percentage (`:percent`) or
in dB (`:dB`). Inputs must be equal-length, nonempty, and finite; the ideal
symbols must have nonzero RMS amplitude.
"""
function evm end

function biterr(tx_bits::AbstractVector{Bool}, rx_bits::AbstractVector{Bool}; mode=:bit, M=nothing)
    bit_count = _validate_metric_vector_lengths(tx_bits, rx_bits, "tx_bits", "rx_bits")
    bits_per_symbol = _validate_biterr_configuration(mode, M, bit_count)
    error_count = 0

    if mode === :bit
        for (tx_bit, rx_bit) in zip(tx_bits, rx_bits)
            error_count += tx_bit != rx_bit
        end
        return error_count, error_count / bit_count
    end

    group_size = 0
    group_has_error = false
    for (tx_bit, rx_bit) in zip(tx_bits, rx_bits)
        group_size += 1
        group_has_error |= tx_bit != rx_bit
        if group_size == bits_per_symbol
            error_count += group_has_error
            group_size = 0
            group_has_error = false
        end
    end
    return error_count, error_count / (bit_count ÷ bits_per_symbol)
end

function evm(ideal_symbols::AbstractVector{<:Complex}, rx_symbols::AbstractVector{<:Complex}; mode=:percent)
    _validate_metric_vector_lengths(ideal_symbols, rx_symbols, "ideal_symbols", "rx_symbols")
    _validate_evm_mode(mode)
    work_type = _evm_work_type(ideal_symbols, rx_symbols)

    if work_type === Float64
        result = _evm_accumulate(ideal_symbols, rx_symbols, Float64)
        if !isnothing(result)
            return _evm_result(result..., mode)
        end

        # Finite Float64 inputs can still overflow in `received - ideal`.
        # Recompute only that exceptional path in BigFloat; ordinary DSP calls
        # keep the Float64, allocation-free path above.
        return _evm_bigfloat_result(ideal_symbols, rx_symbols, mode, Float64)
    end

    return _evm_bigfloat_result(ideal_symbols, rx_symbols, mode, BigFloat)
end

function _evm_accumulate(ideal_symbols::AbstractVector{<:Complex}, rx_symbols::AbstractVector{<:Complex}, ::Type{T}) where {T<:AbstractFloat}
    fast = _evm_accumulate_fast(ideal_symbols, rx_symbols, T)
    isnothing(fast) || return fast
    return _evm_accumulate_scaled(ideal_symbols, rx_symbols, T)
end

function _evm_accumulate_fast(ideal_symbols::AbstractVector{<:Complex}, rx_symbols::AbstractVector{<:Complex}, ::Type{T}) where {T<:AbstractFloat}
    n = length(ideal_symbols)
    ideal_first = firstindex(ideal_symbols)
    rx_first = firstindex(rx_symbols)
    ideal_sumsq = zero(T)
    error_sumsq = zero(T)
    @inbounds @simd for offset in 0:(n - 1)
        ideal = ideal_symbols[ideal_first + offset]
        received = rx_symbols[rx_first + offset]
        ideal_real = T(real(ideal))
        ideal_imag = T(imag(ideal))
        error_real = T(real(received)) - ideal_real
        error_imag = T(imag(received)) - ideal_imag
        ideal_sumsq += ideal_real * ideal_real + ideal_imag * ideal_imag
        error_sumsq += error_real * error_real + error_imag * error_imag
    end

    isfinite(ideal_sumsq) && isfinite(error_sumsq) || return nothing
    ideal_scale = iszero(ideal_sumsq) ? zero(T) : one(T)
    error_scale = iszero(error_sumsq) ? zero(T) : one(T)
    return ideal_scale, ideal_sumsq, error_scale, error_sumsq
end

function _evm_accumulate_scaled(ideal_symbols::AbstractVector{<:Complex}, rx_symbols::AbstractVector{<:Complex}, ::Type{T}) where {T<:AbstractFloat}
    n = length(ideal_symbols)
    ideal_first = firstindex(ideal_symbols)
    rx_first = firstindex(rx_symbols)

    ideal_scale = zero(T)
    error_scale = zero(T)
    @inbounds for offset in 0:(n - 1)
        ideal = ideal_symbols[ideal_first + offset]
        received = rx_symbols[rx_first + offset]
        ideal_real = T(real(ideal))
        ideal_imag = T(imag(ideal))
        error_real = T(real(received)) - ideal_real
        error_imag = T(imag(received)) - ideal_imag
        isfinite(error_real) && isfinite(error_imag) || return nothing
        ideal_scale = max(ideal_scale, hypot(ideal_real, ideal_imag))
        error_scale = max(error_scale, hypot(error_real, error_imag))
    end

    ideal_sumsq = zero(T)
    error_sumsq = zero(T)
    @inbounds @simd for offset in 0:(n - 1)
        ideal = ideal_symbols[ideal_first + offset]
        received = rx_symbols[rx_first + offset]
        ideal_real = T(real(ideal))
        ideal_imag = T(imag(ideal))
        error_real = T(real(received)) - ideal_real
        error_imag = T(imag(received)) - ideal_imag
        ideal_ratio = iszero(ideal_scale) ? zero(T) : hypot(ideal_real, ideal_imag) / ideal_scale
        error_ratio = iszero(error_scale) ? zero(T) : hypot(error_real, error_imag) / error_scale
        ideal_sumsq += ideal_ratio * ideal_ratio
        error_sumsq += error_ratio * error_ratio
    end

    return ideal_scale, ideal_sumsq, error_scale, error_sumsq
end

function _evm_result(ideal_scale::T, ideal_sumsq::T, error_scale::T, error_sumsq::T, mode) where {T<:AbstractFloat}
    iszero(ideal_scale) && throw(ArgumentError("ideal_symbols must have nonzero RMS amplitude"))

    if mode === :dB
        iszero(error_scale) && return T(-Inf)
        return T(20) * (log10(error_scale) - log10(ideal_scale)) + T(10) * (log10(error_sumsq) - log10(ideal_sumsq))
    end

    ratio = (error_scale / ideal_scale) * sqrt(error_sumsq / ideal_sumsq)
    if iszero(ratio) && !iszero(error_scale)
        throw(ArgumentError("EVM percentage is too small to represent in the working precision"))
    end
    isfinite(ratio) || throw(ArgumentError("normalized EVM is not representable in the working precision"))

    percentage = T(100) * ratio
    isfinite(percentage) || throw(ArgumentError("EVM percentage is not representable in the working precision"))
    return percentage
end

function _evm_bigfloat_result(ideal_symbols::AbstractVector{<:Complex}, rx_symbols::AbstractVector{<:Complex}, mode, ::Type{BigFloat})
    result = _evm_accumulate(ideal_symbols, rx_symbols, BigFloat)
    return _evm_result(result..., mode)
end

function _evm_bigfloat_result(ideal_symbols::AbstractVector{<:Complex}, rx_symbols::AbstractVector{<:Complex}, mode, ::Type{Float64})
    result = setprecision(BigFloat, max(precision(BigFloat), 256)) do
        _evm_accumulate(ideal_symbols, rx_symbols, BigFloat)
    end
    big_result = _evm_result(result..., mode)
    float_result = Float64(big_result)
    if mode === :percent
        isfinite(float_result) || throw(ArgumentError("EVM percentage is not representable in the working precision"))
    else
        isfinite(float_result) || throw(ArgumentError("EVM dB is not representable in the working precision"))
    end
    return float_result
end

function _validate_metric_vector_lengths(first::AbstractVector, second::AbstractVector, first_name::AbstractString, second_name::AbstractString)
    length(first) == length(second) || throw(ArgumentError("$first_name and $second_name must have equal lengths"))
    isempty(first) && throw(ArgumentError("$first_name and $second_name must be nonempty"))
    return length(first)
end

function _metric_bits_per_symbol(M)
    M isa Integer || throw(ArgumentError("M must be an integer"))
    M isa Bool && throw(ArgumentError("M must be an integer, not Bool"))
    M >= 2 || throw(ArgumentError("M must be at least 2"))
    M <= typemax(Int) || throw(ArgumentError("M is too large"))

    order = Int(M)
    iszero(order & (order - 1)) || throw(ArgumentError("M must be a power of two"))
    return trailing_zeros(order)
end

function _validate_biterr_configuration(mode, M, bit_length::Integer)
    mode isa Symbol || throw(ArgumentError("mode must be :bit or :symbol"))
    bit_length > 0 || throw(ArgumentError("bit vectors must be nonempty"))

    if mode === :bit
        isnothing(M) || throw(ArgumentError("M is only valid when mode=:symbol"))
        return 1
    elseif mode === :symbol
        isnothing(M) && throw(ArgumentError("M is required when mode=:symbol"))
        bits_per_symbol = _metric_bits_per_symbol(M)
        bit_length % bits_per_symbol == 0 || throw(ArgumentError("bit length must be divisible by log2(M) in symbol mode"))
        return bits_per_symbol
    end

    throw(ArgumentError("mode must be :bit or :symbol"))
end

function _validate_evm_mode(mode)
    mode isa Symbol || throw(ArgumentError("mode must be :percent or :dB"))
    mode in (:percent, :dB) || throw(ArgumentError("mode must be :percent or :dB"))
    return mode
end

function _evm_work_type(ideal_symbols::AbstractVector{<:Complex}, rx_symbols::AbstractVector{<:Complex})
    work_type = Float64
    for (ideal, received) in zip(ideal_symbols, rx_symbols)
        isfinite(ideal) || throw(ArgumentError("ideal_symbols must contain only finite values"))
        isfinite(received) || throw(ArgumentError("rx_symbols must contain only finite values"))
        if !(real(ideal) isa AbstractFloat && imag(ideal) isa AbstractFloat && real(received) isa AbstractFloat && imag(received) isa AbstractFloat) ||
           real(ideal) isa BigFloat || imag(ideal) isa BigFloat || real(received) isa BigFloat || imag(received) isa BigFloat
            work_type = BigFloat
        end
    end
    return work_type
end
