"""
    pskmod(bits, M; gray=true, phase_offset=0.0)

Map MSB-first bits to a unit-modulus M-PSK constellation.
"""
function pskmod end

"""
    pskdemod(symbols, M; gray=true, phase_offset=0.0)

Hard-demodulate unit-modulus M-PSK symbols into MSB-first bits.
"""
function pskdemod end

function pskmod(bits::AbstractVector{<:Integer}, M::Integer; gray::Bool=true, phase_offset::Real=0.0)
    parameters = _psk_parameters(M)
    phase_offset isa Bool && throw(ArgumentError("phase_offset must be a real number, not Bool"))
    isfinite(phase_offset) || throw(ArgumentError("phase_offset must be finite"))
    _validate_bit_vector(bits)
    length(bits) % parameters.bits_per_symbol == 0 || throw(ArgumentError("bit length must be divisible by log2(M)"))

    symbols = Vector{ComplexF64}(undef, length(bits) ÷ parameters.bits_per_symbol)
    constellation = parameters.M <= 65_536 ?
        [cis(phase_offset + 2pi * position / parameters.M) for position in 0:(parameters.M - 1)] : nothing
    input_start = firstindex(bits)
    for symbol_index in eachindex(symbols)
        start = input_start + (symbol_index - 1) * parameters.bits_per_symbol
        label = _bits_to_integer(bits, start, parameters.bits_per_symbol)
        position = gray ? _gray_to_binary(label) : label
        symbols[symbol_index] = isnothing(constellation) ?
            cis(phase_offset + 2pi * position / parameters.M) : constellation[position + 1]
    end
    return symbols
end

function pskdemod(symbols::AbstractVector{<:Complex}, M::Integer; gray::Bool=true, phase_offset::Real=0.0)
    parameters = _psk_parameters(M)
    phase_offset isa Bool && throw(ArgumentError("phase_offset must be a real number, not Bool"))
    isfinite(phase_offset) || throw(ArgumentError("phase_offset must be finite"))
    length(symbols) <= typemax(Int) ÷ parameters.bits_per_symbol || throw(ArgumentError("output bit length is too large"))
    all(isfinite, symbols) || throw(ArgumentError("symbols must contain only finite values"))

    bits = falses(length(symbols) * parameters.bits_per_symbol)
    rotation = cis(-phase_offset)
    rotation8 = cis(-phase_offset - pi / 8)
    rotation16 = cis(-phase_offset - pi / 16)
    for (symbol_index, symbol) in enumerate(symbols)
        position = _psk_demodulate_position(symbol, parameters.M, phase_offset, rotation, rotation8, rotation16)
        label = gray ? _binary_to_gray(position) : position
        start = (symbol_index - 1) * parameters.bits_per_symbol + 1
        _integer_to_bits!(bits, start, parameters.bits_per_symbol, label)
    end
    return bits
end

function _psk_parameters(M::Integer)
    M isa Bool && throw(ArgumentError("M must be an integer, not Bool"))
    M >= 2 || throw(ArgumentError("M must be at least 2"))
    M <= typemax(Int) || throw(ArgumentError("M is too large"))

    order = Int(M)
    iszero(order & (order - 1)) || throw(ArgumentError("M must be a power of two"))

    return (M=order, bits_per_symbol=trailing_zeros(order))
end

function _psk_nearest_position(symbol::Complex, M::Int, phase_offset::Real)
    # Precompute M/(2pi) so exact decision-boundary behavior is stable.
    norm_factor = M / (2pi)
    raw_position = (angle(symbol) - phase_offset) * norm_factor
    return mod(round(Int, raw_position, RoundNearestTiesAway), M)
end

const _PSK_FAST_MARGIN = 1.0e-9

# Only M=2/M=4 have decision boundaries that fall on coordinate axes/diagonals
# after rotating by -phase_offset, letting the position be read off from the
# sign/magnitude comparison of the rotated real/imaginary parts with zero
# trig calls. Returns `nothing` (defer to the exact `_psk_nearest_position`
# atan2 path) whenever the rotated point is close enough to a boundary that
# the two floating-point paths (rotate-then-compare vs. atan2-then-round)
# aren't guaranteed to agree, including the degenerate zero-magnitude case.
function _psk_fast_position(symbol::Complex, M::Int, rotation::Complex)
    rotated = symbol * rotation
    re = real(rotated)
    im = imag(rotated)
    mag2 = re * re + im * im
    if M == 2
        re * re <= _PSK_FAST_MARGIN^2 * mag2 && return nothing
        return re > 0 ? 0 : 1
    else
        diff = re * re - im * im
        abs(diff) <= _PSK_FAST_MARGIN * mag2 && return nothing
        return abs(re) > abs(im) ? (re > 0 ? 0 : 2) : (im > 0 ? 1 : 3)
    end
end


# M=8's true decision boundaries (odd multiples of pi/8) don't fall on the
# axes/diagonals the M=2/M=4 check above uses -- they sit exactly midway
# between them. Rotating by one additional -pi/8 shifts every M=8 boundary
# onto those same axes/diagonals, reducing the problem to the M=4 pattern
# above (now covering all 8 axis/diagonal lines instead of just 4, hence
# three margin checks instead of one) and an extra +1 sector-index shift
# to undo the pi/8 rotation.
function _psk_octant_raw(re8, im8)
    return re8 >= 0 ? (im8 >= 0 ? (re8 >= im8 ? 0 : 1) : (re8 >= -im8 ? 7 : 6)) :
                      (im8 >= 0 ? (-re8 >= im8 ? 3 : 2) : (-re8 >= -im8 ? 4 : 5))
end

function _psk_octant_position(symbol::Complex, rotation8::Complex)
    rotated8 = symbol * rotation8
    re8 = real(rotated8)
    im8 = imag(rotated8)
    mag2 = re8 * re8 + im8 * im8
    re8 * re8 <= _PSK_FAST_MARGIN^2 * mag2 && return nothing
    im8 * im8 <= _PSK_FAST_MARGIN^2 * mag2 && return nothing
    diff = re8 * re8 - im8 * im8
    abs(diff) <= _PSK_FAST_MARGIN * mag2 && return nothing
    return mod(_psk_octant_raw(re8, im8) + 1, 8)
end

const _PSK_OCTANT_ROTATIONS = [cis(-(j * pi) / 4) for j in 0:7]

function _psk_hexadecant_position(symbol::Complex, rotation16::Complex)
    rotated16 = symbol * rotation16
    re16 = real(rotated16)
    im16 = imag(rotated16)
    mag2 = re16 * re16 + im16 * im16
    re16 * re16 <= _PSK_FAST_MARGIN^2 * mag2 && return nothing
    im16 * im16 <= _PSK_FAST_MARGIN^2 * mag2 && return nothing
    diff = re16 * re16 - im16 * im16
    abs(diff) <= _PSK_FAST_MARGIN * mag2 && return nothing
    j = _psk_octant_raw(re16, im16)
    w = rotated16 * (@inbounds _PSK_OCTANT_ROTATIONS[j + 1])
    bisector = imag(w) - (sqrt(2.0) - 1.0) * real(w)
    bisector * bisector <= _PSK_FAST_MARGIN^2 * mag2 && return nothing
    return mod(2 * j + (bisector > 0 ? 1 : 0) + 1, 16)
end

function _psk_demodulate_position(symbol::Complex, M::Int, phase_offset::Real, rotation::Complex, rotation8::Complex, rotation16::Complex)
    if M == 2 || M == 4
        fast = _psk_fast_position(symbol, M, rotation)
        isnothing(fast) || return fast
    elseif M == 8
        fast = _psk_octant_position(symbol, rotation8)
        isnothing(fast) || return fast
    elseif M == 16
        fast = _psk_hexadecant_position(symbol, rotation16)
        isnothing(fast) || return fast
    end
    return _psk_nearest_position(symbol, M, phase_offset)
end
