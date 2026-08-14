const _SOFT_DEMOD_MAX_ORDER = 4096

"""
    qamsoftdemod(symbols, M; gray=true, noise_variance=1.0)

Return max-log bit log-likelihood ratios for a unit-average-power QAM
constellation. Results are MSB-first and positive values favor bit `0`.
`noise_variance` is the complex AWGN variance `E[abs2(noise)]` at the
demodulator input.
"""
function qamsoftdemod end

"""
    psksoftdemod(symbols, M; gray=true, phase_offset=0.0, noise_variance=1.0)

Return max-log bit log-likelihood ratios for a unit-modulus M-PSK
constellation. Results are MSB-first and positive values favor bit `0`.
`noise_variance` is the complex AWGN variance `E[abs2(noise)]` at the
demodulator input.
"""
function psksoftdemod end

function qamsoftdemod(symbols::AbstractVector{<:Complex}, M::Integer;
                      gray::Bool=true, noise_variance::Real=1.0)
    parameters = _qam_parameters(M)
    variance = _soft_demod_variance(noise_variance)
    _validate_soft_demod_symbols(symbols, parameters.bits_per_symbol)
    isempty(symbols) && return Float64[]
    parameters.M <= _SOFT_DEMOD_MAX_ORDER ||
        throw(ArgumentError("soft demodulation supports M <= $_SOFT_DEMOD_MAX_ORDER"))

    parameters.kind === :square &&
        return _qam_square_maxlog_llrs(symbols, parameters, gray, variance)
    labels = _soft_demod_label_bits(parameters.M, parameters.bits_per_symbol)
    constellation = qammod(labels, parameters.M; gray)
    return _maxlog_llrs(symbols, constellation, parameters.bits_per_symbol, variance)
end

function psksoftdemod(symbols::AbstractVector{<:Complex}, M::Integer;
                      gray::Bool=true, phase_offset::Real=0.0, noise_variance::Real=1.0)
    parameters = _psk_parameters(M)
    phase_offset isa Bool && throw(ArgumentError("phase_offset must be a real number, not Bool"))
    isfinite(phase_offset) || throw(ArgumentError("phase_offset must be finite"))
    variance = _soft_demod_variance(noise_variance)
    _validate_soft_demod_symbols(symbols, parameters.bits_per_symbol)
    isempty(symbols) && return Float64[]
    parameters.M <= _SOFT_DEMOD_MAX_ORDER ||
        throw(ArgumentError("soft demodulation supports M <= $_SOFT_DEMOD_MAX_ORDER"))

    labels = _soft_demod_label_bits(parameters.M, parameters.bits_per_symbol)
    constellation = pskmod(labels, parameters.M; gray, phase_offset)
    return _maxlog_llrs(symbols, constellation, parameters.bits_per_symbol, variance, true)
end

function _soft_demod_variance(noise_variance::Real)
    noise_variance isa Bool && throw(ArgumentError("noise_variance must be a positive real number, not Bool"))
    variance = Float64(noise_variance)
    (isfinite(variance) && variance > 0) ||
        throw(ArgumentError("noise_variance must be finite, positive, and representable as Float64"))
    return variance
end

function _validate_soft_demod_symbols(symbols, bits_per_symbol::Int)
    length(symbols) <= typemax(Int) ÷ bits_per_symbol ||
        throw(ArgumentError("output LLR length is too large"))
    all(isfinite, symbols) || throw(ArgumentError("symbols must contain only finite values"))
    return nothing
end

function _soft_demod_label_bits(M::Int, bits_per_symbol::Int)
    labels = falses(M * bits_per_symbol)
    for label in 0:(M - 1)
        _integer_to_bits!(labels, label * bits_per_symbol + 1, bits_per_symbol, label)
    end
    return labels
end

function _qam_square_maxlog_llrs(symbols, parameters, gray::Bool,
                                  noise_variance::Float64)
    bits_per_symbol = parameters.bits_per_symbol
    axis_bits = parameters.axis_bits
    llrs = Vector{Float64}(undef, length(symbols) * bits_per_symbol)
    minimum_zero = Vector{Float64}(undef, axis_bits)
    minimum_one = Vector{Float64}(undef, axis_bits)

    for (symbol_index, symbol) in enumerate(symbols)
        start = (symbol_index - 1) * bits_per_symbol + 1
        real_value = Float64(real(symbol))
        imag_value = Float64(imag(symbol))
        if !isfinite(real_value) || !isfinite(imag_value) ||
           max(abs(real_value), abs(imag_value)) > sqrt(floatmax(Float64))
            _qam_square_symbol_llrs_big!(llrs, start, symbol, parameters, gray,
                                         noise_variance)
            continue
        end

        _qam_square_axis_llrs!(llrs, start, real_value, parameters, gray, false,
                               noise_variance, minimum_zero, minimum_one)
        _qam_square_axis_llrs!(llrs, start + axis_bits, imag_value, parameters,
                               gray, true, noise_variance, minimum_zero, minimum_one)
    end
    return llrs
end

function _qam_square_axis_llrs!(llrs, start::Int, value::Float64, parameters,
                                 gray::Bool, reverse::Bool, noise_variance::Float64,
                                 minimum_zero, minimum_one)
    axis_bits = parameters.axis_bits
    side_length = parameters.side_length
    fill!(minimum_zero, Inf)
    fill!(minimum_one, Inf)

    for grid_position in 0:(side_length - 1)
        level = parameters.scale * _qam_axis_level(grid_position, side_length)
        score = level * level - 2 * value * level
        label_position = reverse ? side_length - 1 - grid_position : grid_position
        label = gray ? _binary_to_gray(label_position) : label_position
        for bit_index in 1:axis_bits
            bit = (label >> (axis_bits - bit_index)) & 1
            if bit == 0
                minimum_zero[bit_index] = min(minimum_zero[bit_index], score)
            else
                minimum_one[bit_index] = min(minimum_one[bit_index], score)
            end
        end
    end

    for bit_index in 1:axis_bits
        llrs[start + bit_index - 1] =
            (minimum_one[bit_index] - minimum_zero[bit_index]) / noise_variance
    end
    return llrs
end

function _qam_square_symbol_llrs_big!(llrs, start::Int, symbol, parameters,
                                       gray::Bool, noise_variance::Float64)
    input_precision = real(symbol) isa BigFloat ? precision(real(symbol)) : 0
    return setprecision(BigFloat, max(2048, input_precision + 64)) do
        axis_bits = parameters.axis_bits
        side_length = parameters.side_length
        scale = BigFloat(parameters.scale)
        variance = BigFloat(noise_variance)
        limit = BigFloat(floatmax(Float64))
        values = (BigFloat(real(symbol)), BigFloat(imag(symbol)))

        for axis in 1:2
            minimum_zero = fill(BigFloat(Inf), axis_bits)
            minimum_one = fill(BigFloat(Inf), axis_bits)
            reverse = axis == 2
            for grid_position in 0:(side_length - 1)
                level = scale * _qam_axis_level(grid_position, side_length)
                score = level^2 - 2 * values[axis] * level
                label_position = reverse ? side_length - 1 - grid_position : grid_position
                label = gray ? _binary_to_gray(label_position) : label_position
                for bit_index in 1:axis_bits
                    bit = (label >> (axis_bits - bit_index)) & 1
                    if bit == 0
                        minimum_zero[bit_index] = min(minimum_zero[bit_index], score)
                    else
                        minimum_one[bit_index] = min(minimum_one[bit_index], score)
                    end
                end
            end
            output_start = start + (axis - 1) * axis_bits
            for bit_index in 1:axis_bits
                llr = (minimum_one[bit_index] - minimum_zero[bit_index]) / variance
                llrs[output_start + bit_index - 1] = Float64(clamp(llr, -limit, limit))
            end
        end
        return llrs
    end
end

function _maxlog_llrs(symbols, constellation::Vector{ComplexF64}, bits_per_symbol::Int,
                       noise_variance::Float64, constant_energy::Bool=false)
    llrs = Vector{Float64}(undef, length(symbols) * bits_per_symbol)
    minimum_zero = Vector{Float64}(undef, bits_per_symbol)
    minimum_one = Vector{Float64}(undef, bits_per_symbol)
    for (symbol_index, symbol) in enumerate(symbols)
        start = (symbol_index - 1) * bits_per_symbol + 1
        _maxlog_symbol_llrs!(llrs, start, symbol, constellation, bits_per_symbol,
                             noise_variance, minimum_zero, minimum_one, constant_energy)
    end
    return llrs
end

function _maxlog_symbol_llrs!(llrs::Vector{Float64}, start::Int, symbol,
                               constellation::Vector{ComplexF64}, bits_per_symbol::Int,
                               noise_variance::Float64, minimum_zero::Vector{Float64},
                               minimum_one::Vector{Float64}, constant_energy::Bool=false)
    real_value = Float64(real(symbol))
    imag_value = Float64(imag(symbol))
    if !isfinite(real_value) || !isfinite(imag_value) ||
       max(abs(real_value), abs(imag_value)) > sqrt(floatmax(Float64))
        return _maxlog_symbol_llrs_big!(llrs, start, symbol, constellation, bits_per_symbol,
                                        noise_variance, constant_energy)
    end

    fill!(minimum_zero, Inf)
    fill!(minimum_one, Inf)
    for (index, point) in enumerate(constellation)
        score = (constant_energy ? 0.0 : abs2(point)) -
                2 * (real_value * real(point) + imag_value * imag(point))
        isfinite(score) ||
            return _maxlog_symbol_llrs_big!(llrs, start, symbol, constellation, bits_per_symbol,
                                            noise_variance, constant_energy)
        label = index - 1
        for bit_index in 1:bits_per_symbol
            bit = (label >> (bits_per_symbol - bit_index)) & 1
            if bit == 0
                minimum_zero[bit_index] = min(minimum_zero[bit_index], score)
            else
                minimum_one[bit_index] = min(minimum_one[bit_index], score)
            end
        end
    end

    for bit_index in 1:bits_per_symbol
        llr = (minimum_one[bit_index] - minimum_zero[bit_index]) / noise_variance
        isfinite(llr) ||
            return _maxlog_symbol_llrs_big!(llrs, start, symbol, constellation, bits_per_symbol,
                                            noise_variance, constant_energy)
        llrs[start + bit_index - 1] = llr
    end
    return llrs
end

function _maxlog_symbol_llrs_big!(llrs::Vector{Float64}, start::Int, symbol,
                                   constellation::Vector{ComplexF64}, bits_per_symbol::Int,
                                   noise_variance::Float64, constant_energy::Bool=false)
    input_precision = real(symbol) isa BigFloat ? precision(real(symbol)) : 0
    return setprecision(BigFloat, max(2048, input_precision + 64)) do
        real_value = BigFloat(real(symbol))
        imag_value = BigFloat(imag(symbol))
        minimum_zero = fill(BigFloat(Inf), bits_per_symbol)
        minimum_one = fill(BigFloat(Inf), bits_per_symbol)

        for (index, point) in enumerate(constellation)
            point_real = BigFloat(real(point))
            point_imag = BigFloat(imag(point))
            score = (constant_energy ? 0 : point_real^2 + point_imag^2) -
                    2 * (real_value * point_real + imag_value * point_imag)
            label = index - 1
            for bit_index in 1:bits_per_symbol
                bit = (label >> (bits_per_symbol - bit_index)) & 1
                if bit == 0
                    minimum_zero[bit_index] = min(minimum_zero[bit_index], score)
                else
                    minimum_one[bit_index] = min(minimum_one[bit_index], score)
                end
            end
        end

        limit = BigFloat(floatmax(Float64))
        variance = BigFloat(noise_variance)
        for bit_index in 1:bits_per_symbol
            llr = (minimum_one[bit_index] - minimum_zero[bit_index]) / variance
            llrs[start + bit_index - 1] = Float64(clamp(llr, -limit, limit))
        end
        return llrs
    end
end
