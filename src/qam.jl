"""
    qammod(bits, M; gray=true)

Map MSB-first bits to a unit-average-power square, rectangular, or cross-QAM
constellation. Supported orders are 4, 8, 16, 32, 64, 128, 256, and 512.
"""
function qammod end

"""
    qamdemod(symbols, M; gray=true)

Hard-demodulate a supported unit-average-power QAM constellation into
MSB-first bits.
"""
function qamdemod end

function qammod(bits::AbstractVector{<:Integer}, M::Integer; gray::Bool=true)
    parameters = _qam_parameters(M)
    _validate_bit_vector(bits)
    length(bits) % parameters.bits_per_symbol == 0 || throw(ArgumentError("bit length must be divisible by log2(M)"))

    symbols = Vector{ComplexF64}(undef, length(bits) ÷ parameters.bits_per_symbol)
    _qam_modulate_dispatch!(symbols, bits, parameters, gray, Val(parameters.kind))
    return symbols
end

function qamdemod(symbols::AbstractVector{<:Complex}, M::Integer; gray::Bool=true)
    parameters = _qam_parameters(M)
    length(symbols) <= typemax(Int) ÷ parameters.bits_per_symbol || throw(ArgumentError("output bit length is too large"))
    all(isfinite, symbols) || throw(ArgumentError("symbols must contain only finite values"))

    bits = falses(length(symbols) * parameters.bits_per_symbol)
    _qam_demodulate_dispatch!(bits, symbols, parameters, gray, Val(parameters.kind))
    return bits
end

function _qam_modulate_dispatch!(symbols, bits, parameters, gray::Bool, ::Val{:square})
    input_start = firstindex(bits)
    for symbol_index in eachindex(symbols)
        start = input_start + (symbol_index - 1) * parameters.bits_per_symbol
        i_label = _bits_to_integer(bits, start, parameters.axis_bits)
        q_label = _bits_to_integer(bits, start + parameters.axis_bits, parameters.axis_bits)
        i_position = gray ? _gray_to_binary(i_label) : i_label
        q_position = gray ? _gray_to_binary(q_label) : q_label
        i_level = _qam_axis_level(i_position, parameters.side_length)
        q_level = _qam_axis_level(parameters.side_length - 1 - q_position, parameters.side_length)
        symbols[symbol_index] = parameters.scale * complex(i_level, q_level)
    end
end

function _qam_modulate_dispatch!(symbols, bits, parameters, gray::Bool, ::Val{:nonsquare})
    input_start = firstindex(bits)
    table = _qam_nonsquare_table(parameters.M, gray)
    for symbol_index in eachindex(symbols)
        start = input_start + (symbol_index - 1) * parameters.bits_per_symbol
        label = _bits_to_integer(bits, start, parameters.bits_per_symbol)
        i_level, q_level = table[label + 1]
        symbols[symbol_index] = parameters.scale * complex(i_level, q_level)
    end
end

function _qam_demodulate_dispatch!(bits, symbols, parameters, gray::Bool, ::Val{:square})
    # Keep the normalization order explicit for stable boundary decisions.
    inverse_scale = sqrt(2.0 * (Float64(parameters.M) - 1.0) / 3.0)
    for (symbol_index, symbol) in enumerate(symbols)
        i_position = _qam_nearest_axis_position(real(symbol) * inverse_scale, parameters.side_length)
        q_grid_position = _qam_nearest_axis_position(imag(symbol) * inverse_scale, parameters.side_length)
        i_label = gray ? _binary_to_gray(i_position) : i_position
        q_position = parameters.side_length - 1 - q_grid_position
        q_label = gray ? _binary_to_gray(q_position) : q_position
        start = (symbol_index - 1) * parameters.bits_per_symbol + 1
        _integer_to_bits!(bits, start, parameters.axis_bits, i_label)
        _integer_to_bits!(bits, start + parameters.axis_bits, parameters.axis_bits, q_label)
    end
end

function _qam_demodulate_dispatch!(bits, symbols, parameters, gray::Bool, ::Val{:nonsquare})
    # Descale by sqrt(averagePower), then choose the minimum-distance point.
    # Ties use the first bin-order position, independent of Gray labeling.
    # Most points are handled in
    # O(1) by `_qam_nonsquare_fast_bin_label` (provably correct away from
    # the removed corners and from genuine ties); anything it can't decide
    # falls back to the exact O(M) brute-force search over the bin table.
    inverse_scale = sqrt(parameters.average_power)
    bin_table = _qam_nonsquare_table(parameters.M, false)
    bin_to_gray = _qam_nonsquare_bin_to_gray_table(parameters.M)
    for (symbol_index, symbol) in enumerate(symbols)
        i_value = real(symbol) * inverse_scale
        q_value = imag(symbol) * inverse_scale
        if !isfinite(i_value) || !isfinite(q_value)
            bin_position = _qam_nonsquare_nearest_symbol_label(real(symbol), imag(symbol), parameters.scale, bin_table)
        else
            fast_label = _qam_nonsquare_fast_bin_label(i_value, q_value, parameters)
            bin_position = isnothing(fast_label) ? _qam_nonsquare_nearest_label(i_value, q_value, bin_table) : fast_label
        end
        label = gray ? bin_to_gray[bin_position + 1] : bin_position
        start = (symbol_index - 1) * parameters.bits_per_symbol + 1
        _integer_to_bits!(bits, start, parameters.bits_per_symbol, label)
    end
end

function _qam_nonsquare_nearest_symbol_label(i_value::Real, q_value::Real, scale::Float64,
                                              table::Vector{Tuple{Int,Int}})
    i_precision = i_value isa BigFloat ? precision(i_value) : 0
    q_precision = q_value isa BigFloat ? precision(q_value) : 0
    return setprecision(BigFloat, max(2048, i_precision + 64, q_precision + 64)) do
        i_big = BigFloat(i_value)
        q_big = BigFloat(q_value)
        scale_big = BigFloat(scale)
        best_label = 0
        best_distance = Inf
        for (index, (i_level, q_level)) in enumerate(table)
            distance = abs2(i_big - scale_big * i_level) + abs2(q_big - scale_big * q_level)
            if distance < best_distance
                best_distance = distance
                best_label = index - 1
            end
        end
        return best_label
    end
end

function _qam_parameters(M::Integer)
    M isa Bool && throw(ArgumentError("M must be an integer, not Bool"))
    M >= 4 || throw(ArgumentError("M must be at least 4"))
    M <= typemax(Int) || throw(ArgumentError("M is too large"))

    order = Int(M)
    iszero(order & (order - 1)) || throw(ArgumentError("M must be a power of two"))

    bits_per_symbol = trailing_zeros(order)
    if iseven(bits_per_symbol)
        return _square_qam_parameters(order, bits_per_symbol)
    else
        order in (8, 32, 128, 512) || throw(ArgumentError("M=$order is not a supported QAM order (square: 4,16,64,256,...; cross/rectangular: 8,32,128,512)"))
        return _nonsquare_qam_parameters(order, bits_per_symbol)
    end
end

function _square_qam_parameters(order::Int, bits_per_symbol::Int)
    axis_bits = bits_per_symbol ÷ 2
    side_length = 1 << axis_bits
    return (
        M=order,
        bits_per_symbol=bits_per_symbol,
        kind=:square,
        axis_bits=axis_bits,
        side_length=side_length,
        scale=sqrt(3.0 / (2.0 * (Float64(order) - 1.0))),
    )
end

function _nonsquare_qam_parameters(order::Int, bits_per_symbol::Int)
    # Closed-form average power for minimum-distance-2 constellations. M=8 is
    # rectangular; larger supported non-square orders are corner-cut crosses.
    average_power = order == 8 ? (5.0 * order / 4.0 - 1.0) * 2.0 / 3.0 : (31.0 * order / 32.0 - 1.0) * 2.0 / 3.0
    if order == 8
        i_side, q_side, corner_size = 4, 2, 0
    else
        half_side = 1 << ((bits_per_symbol - 1) ÷ 2)
        grid_side = (3 * half_side) ÷ 2
        i_side, q_side, corner_size = grid_side, grid_side, half_side ÷ 4
    end
    return (
        M=order,
        bits_per_symbol=bits_per_symbol,
        kind=:nonsquare,
        i_side=i_side,
        q_side=q_side,
        corner_size=corner_size,
        average_power=average_power,
        scale=sqrt(1.0 / average_power),
    )
end

_qam_axis_level(position::Int, side_length::Int) = 2 * position - side_length + 1

function _qam_nearest_axis_position(value::Real, side_length::Int)
    position = (value + side_length - 1) / 2
    position <= 0 && return 0
    position >= side_length - 1 && return side_length - 1
    return floor(Int, position + 1 / 2)
end

@inline function _qam_nonsquare_table(M::Int, gray::Bool)
    if M == 8
        return gray ? _QAM_RECTANGULAR_8_GRAY : _QAM_RECTANGULAR_8_BIN
    elseif M == 32
        return gray ? _QAM_CROSS_32_GRAY : _QAM_CROSS_32_BIN
    elseif M == 128
        return gray ? _QAM_CROSS_128_GRAY : _QAM_CROSS_128_BIN
    else
        return gray ? _QAM_CROSS_512_GRAY : _QAM_CROSS_512_BIN
    end
end

function _qam_nonsquare_nearest_label(i_value::Real, q_value::Real, table::Vector{Tuple{Int,Int}})
    best_label = 0
    best_distance = Inf
    for (index, (i_level, q_level)) in enumerate(table)
        distance = abs2(i_value - i_level) + abs2(q_value - q_level)
        isfinite(distance) || return _qam_nonsquare_nearest_label_bigfloat(i_value, q_value, table)
        if distance < best_distance
            best_distance = distance
            best_label = index - 1
        end
    end
    return best_label
end

function _qam_nonsquare_nearest_label_bigfloat(i_value::Real, q_value::Real, table::Vector{Tuple{Int,Int}})
    i_precision = i_value isa BigFloat ? precision(i_value) : 0
    q_precision = q_value isa BigFloat ? precision(q_value) : 0
    return setprecision(BigFloat, max(2048, i_precision + 64, q_precision + 64)) do
        i_big = BigFloat(i_value)
        q_big = BigFloat(q_value)
        best_label = 0
        best_distance = Inf
        for (index, (i_level, q_level)) in enumerate(table)
            distance = abs2(i_big - i_level) + abs2(q_big - q_level)
            if distance < best_distance
                best_distance = distance
                best_label = index - 1
            end
        end
        return best_label
    end
end

@inline function _qam_nonsquare_bin_to_gray_table(M::Int)
    if M == 8
        return _QAM_RECTANGULAR_8_BIN_TO_GRAY
    elseif M == 32
        return _QAM_CROSS_32_BIN_TO_GRAY
    elseif M == 128
        return _QAM_CROSS_128_BIN_TO_GRAY
    else
        return _QAM_CROSS_512_BIN_TO_GRAY
    end
end

function _qam_nonsquare_build_bin_to_gray(bin_table::Vector{Tuple{Int,Int}}, gray_table::Vector{Tuple{Int,Int}})
    mapping = Vector{Int}(undef, length(bin_table))
    for (bin_position, level) in enumerate(bin_table)
        mapping[bin_position] = findfirst(==(level), gray_table) - 1
    end
    return mapping
end

function _qam_nonsquare_is_corner(i_position::Int, q_position::Int, i_side::Int, q_side::Int, corner_size::Int)
    i_edge = i_position < corner_size || i_position >= i_side - corner_size
    q_edge = q_position < corner_size || q_position >= q_side - corner_size
    return i_edge && q_edge
end

# The per-axis nearest position independently minimizes each axis' squared
# distance over the *full* uncut grid, so it's a lower bound on total
# distance to any point, real or not. If it isn't a removed corner
# position, it's an actual point achieving that lower bound, hence the true
# global nearest neighbor -- no search needed. Since each axis' distance is
# a strictly increasing parabola away from its own vertex, only the
# immediately adjacent position on each axis can ever tie it, so checking
# the up-to-8 neighboring cells (using the exact same distance formula as
# the ground-truth brute-force search) is a complete tie check, not a
# heuristic one. Falls back to `nothing` (meaning: let the caller run the
# existing exact search) for corners and for any tie.
function _qam_nonsquare_fast_bin_label(i_value::Real, q_value::Real, parameters)
    i_side = parameters.i_side
    q_side = parameters.q_side
    corner_size = parameters.corner_size
    i_position = _qam_nearest_axis_position(i_value, i_side)
    q_position = _qam_nearest_axis_position(q_value, q_side)
    _qam_nonsquare_is_corner(i_position, q_position, i_side, q_side, corner_size) && return nothing

    best_distance = abs2(i_value - _qam_axis_level(i_position, i_side)) + abs2(q_value - _qam_axis_level(q_position, q_side))
    for i_offset in -1:1, q_offset in -1:1
        (i_offset == 0 && q_offset == 0) && continue
        neighbor_i = i_position + i_offset
        neighbor_q = q_position + q_offset
        (0 <= neighbor_i < i_side && 0 <= neighbor_q < q_side) || continue
        _qam_nonsquare_is_corner(neighbor_i, neighbor_q, i_side, q_side, corner_size) && continue
        neighbor_distance = abs2(i_value - _qam_axis_level(neighbor_i, i_side)) + abs2(q_value - _qam_axis_level(neighbor_q, q_side))
        neighbor_distance <= best_distance && return nothing
    end

    index = _qam_nonsquare_position_index(parameters.M)
    return index[i_position * q_side + q_position + 1]
end

function _qam_nonsquare_build_position_index(bin_table::Vector{Tuple{Int,Int}}, i_side::Int, q_side::Int)
    index = fill(-1, i_side * q_side)
    for (bin_label, (i_level, q_level)) in enumerate(bin_table)
        i_position = (i_level + i_side - 1) ÷ 2
        q_position = (q_level + q_side - 1) ÷ 2
        index[i_position * q_side + q_position + 1] = bin_label - 1
    end
    return index
end

@inline function _qam_nonsquare_position_index(M::Int)
    if M == 8
        return _QAM_RECTANGULAR_8_POSITION_INDEX
    elseif M == 32
        return _QAM_CROSS_32_POSITION_INDEX
    elseif M == 128
        return _QAM_CROSS_128_POSITION_INDEX
    else
        return _QAM_CROSS_512_POSITION_INDEX
    end
end

# Label-to-level tables for non-square QAM. M=8 is a 4x2 rectangle; M=32,
# 128, and 512 are corner-cut crosses whose Gray labeling is stored directly.
const _QAM_RECTANGULAR_8_GRAY = [(-3,1),(-3,-1),(-1,1),(-1,-1),(3,1),(3,-1),(1,1),(1,-1)]
const _QAM_RECTANGULAR_8_BIN = [(-3,1),(-3,-1),(-1,1),(-1,-1),(1,1),(1,-1),(3,1),(3,-1)]
const _QAM_CROSS_32_GRAY = [(-3,5),(-1,5),(-3,-5),(-1,-5),(-5,3),(-5,1),(-5,-3),(-5,-1),(-1,3),(-1,1),(-1,-3),(-1,-1),(-3,3),(-3,1),(-3,-3),(-3,-1),(3,5),(1,5),(3,-5),(1,-5),(5,3),(5,1),(5,-3),(5,-1),(1,3),(1,1),(1,-3),(1,-1),(3,3),(3,1),(3,-3),(3,-1)]
const _QAM_CROSS_32_BIN = [(-3,5),(-1,5),(-1,-5),(-3,-5),(-5,3),(-5,1),(-5,-1),(-5,-3),(-3,3),(-3,1),(-3,-1),(-3,-3),(-1,3),(-1,1),(-1,-1),(-1,-3),(1,3),(1,1),(1,-1),(1,-3),(3,3),(3,1),(3,-1),(3,-3),(5,3),(5,1),(5,-1),(5,-3),(3,5),(1,5),(1,-5),(3,-5)]
const _QAM_CROSS_128_GRAY = [(-7,9),(-7,11),(-1,9),(-1,11),(-7,-9),(-7,-11),(-1,-9),(-1,-11),(-5,9),(-5,11),(-3,9),(-3,11),(-5,-9),(-5,-11),(-3,-9),(-3,-11),(-9,7),(-9,5),(-9,1),(-9,3),(-9,-7),(-9,-5),(-9,-1),(-9,-3),(-11,7),(-11,5),(-11,1),(-11,3),(-11,-7),(-11,-5),(-11,-1),(-11,-3),(-1,7),(-1,5),(-1,1),(-1,3),(-1,-7),(-1,-5),(-1,-1),(-1,-3),(-3,7),(-3,5),(-3,1),(-3,3),(-3,-7),(-3,-5),(-3,-1),(-3,-3),(-7,7),(-7,5),(-7,1),(-7,3),(-7,-7),(-7,-5),(-7,-1),(-7,-3),(-5,7),(-5,5),(-5,1),(-5,3),(-5,-7),(-5,-5),(-5,-1),(-5,-3),(7,9),(7,11),(1,9),(1,11),(7,-9),(7,-11),(1,-9),(1,-11),(5,9),(5,11),(3,9),(3,11),(5,-9),(5,-11),(3,-9),(3,-11),(9,7),(9,5),(9,1),(9,3),(9,-7),(9,-5),(9,-1),(9,-3),(11,7),(11,5),(11,1),(11,3),(11,-7),(11,-5),(11,-1),(11,-3),(1,7),(1,5),(1,1),(1,3),(1,-7),(1,-5),(1,-1),(1,-3),(3,7),(3,5),(3,1),(3,3),(3,-7),(3,-5),(3,-1),(3,-3),(7,7),(7,5),(7,1),(7,3),(7,-7),(7,-5),(7,-1),(7,-3),(5,7),(5,5),(5,1),(5,3),(5,-7),(5,-5),(5,-1),(5,-3)]
const _QAM_CROSS_128_BIN = [(-7,9),(-7,11),(-1,11),(-1,9),(-1,-9),(-1,-11),(-7,-11),(-7,-9),(-5,9),(-5,11),(-3,11),(-3,9),(-3,-9),(-3,-11),(-5,-11),(-5,-9),(-11,7),(-11,5),(-11,3),(-11,1),(-11,-1),(-11,-3),(-11,-5),(-11,-7),(-9,7),(-9,5),(-9,3),(-9,1),(-9,-1),(-9,-3),(-9,-5),(-9,-7),(-7,7),(-7,5),(-7,3),(-7,1),(-7,-1),(-7,-3),(-7,-5),(-7,-7),(-5,7),(-5,5),(-5,3),(-5,1),(-5,-1),(-5,-3),(-5,-5),(-5,-7),(-3,7),(-3,5),(-3,3),(-3,1),(-3,-1),(-3,-3),(-3,-5),(-3,-7),(-1,7),(-1,5),(-1,3),(-1,1),(-1,-1),(-1,-3),(-1,-5),(-1,-7),(1,7),(1,5),(1,3),(1,1),(1,-1),(1,-3),(1,-5),(1,-7),(3,7),(3,5),(3,3),(3,1),(3,-1),(3,-3),(3,-5),(3,-7),(5,7),(5,5),(5,3),(5,1),(5,-1),(5,-3),(5,-5),(5,-7),(7,7),(7,5),(7,3),(7,1),(7,-1),(7,-3),(7,-5),(7,-7),(9,7),(9,5),(9,3),(9,1),(9,-1),(9,-3),(9,-5),(9,-7),(11,7),(11,5),(11,3),(11,1),(11,-1),(11,-3),(11,-5),(11,-7),(5,9),(5,11),(3,11),(3,9),(3,-9),(3,-11),(5,-11),(5,-9),(7,9),(7,11),(1,11),(1,9),(1,-9),(1,-11),(7,-11),(7,-9)]
const _QAM_CROSS_512_GRAY = [(-15,17),(-15,19),(-15,23),(-15,21),(-1,17),(-1,19),(-1,23),(-1,21),(-15,-17),(-15,-19),(-15,-23),(-15,-21),(-1,-17),(-1,-19),(-1,-23),(-1,-21),(-13,17),(-13,19),(-13,23),(-13,21),(-3,17),(-3,19),(-3,23),(-3,21),(-13,-17),(-13,-19),(-13,-23),(-13,-21),(-3,-17),(-3,-19),(-3,-23),(-3,-21),(-9,17),(-9,19),(-9,23),(-9,21),(-7,17),(-7,19),(-7,23),(-7,21),(-9,-17),(-9,-19),(-9,-23),(-9,-21),(-7,-17),(-7,-19),(-7,-23),(-7,-21),(-11,17),(-11,19),(-11,23),(-11,21),(-5,17),(-5,19),(-5,23),(-5,21),(-11,-17),(-11,-19),(-11,-23),(-11,-21),(-5,-17),(-5,-19),(-5,-23),(-5,-21),(-17,15),(-17,13),(-17,9),(-17,11),(-17,1),(-17,3),(-17,7),(-17,5),(-17,-15),(-17,-13),(-17,-9),(-17,-11),(-17,-1),(-17,-3),(-17,-7),(-17,-5),(-19,15),(-19,13),(-19,9),(-19,11),(-19,1),(-19,3),(-19,7),(-19,5),(-19,-15),(-19,-13),(-19,-9),(-19,-11),(-19,-1),(-19,-3),(-19,-7),(-19,-5),(-23,15),(-23,13),(-23,9),(-23,11),(-23,1),(-23,3),(-23,7),(-23,5),(-23,-15),(-23,-13),(-23,-9),(-23,-11),(-23,-1),(-23,-3),(-23,-7),(-23,-5),(-21,15),(-21,13),(-21,9),(-21,11),(-21,1),(-21,3),(-21,7),(-21,5),(-21,-15),(-21,-13),(-21,-9),(-21,-11),(-21,-1),(-21,-3),(-21,-7),(-21,-5),(-1,15),(-1,13),(-1,9),(-1,11),(-1,1),(-1,3),(-1,7),(-1,5),(-1,-15),(-1,-13),(-1,-9),(-1,-11),(-1,-1),(-1,-3),(-1,-7),(-1,-5),(-3,15),(-3,13),(-3,9),(-3,11),(-3,1),(-3,3),(-3,7),(-3,5),(-3,-15),(-3,-13),(-3,-9),(-3,-11),(-3,-1),(-3,-3),(-3,-7),(-3,-5),(-7,15),(-7,13),(-7,9),(-7,11),(-7,1),(-7,3),(-7,7),(-7,5),(-7,-15),(-7,-13),(-7,-9),(-7,-11),(-7,-1),(-7,-3),(-7,-7),(-7,-5),(-5,15),(-5,13),(-5,9),(-5,11),(-5,1),(-5,3),(-5,7),(-5,5),(-5,-15),(-5,-13),(-5,-9),(-5,-11),(-5,-1),(-5,-3),(-5,-7),(-5,-5),(-15,15),(-15,13),(-15,9),(-15,11),(-15,1),(-15,3),(-15,7),(-15,5),(-15,-15),(-15,-13),(-15,-9),(-15,-11),(-15,-1),(-15,-3),(-15,-7),(-15,-5),(-13,15),(-13,13),(-13,9),(-13,11),(-13,1),(-13,3),(-13,7),(-13,5),(-13,-15),(-13,-13),(-13,-9),(-13,-11),(-13,-1),(-13,-3),(-13,-7),(-13,-5),(-9,15),(-9,13),(-9,9),(-9,11),(-9,1),(-9,3),(-9,7),(-9,5),(-9,-15),(-9,-13),(-9,-9),(-9,-11),(-9,-1),(-9,-3),(-9,-7),(-9,-5),(-11,15),(-11,13),(-11,9),(-11,11),(-11,1),(-11,3),(-11,7),(-11,5),(-11,-15),(-11,-13),(-11,-9),(-11,-11),(-11,-1),(-11,-3),(-11,-7),(-11,-5),(15,17),(15,19),(15,23),(15,21),(1,17),(1,19),(1,23),(1,21),(15,-17),(15,-19),(15,-23),(15,-21),(1,-17),(1,-19),(1,-23),(1,-21),(13,17),(13,19),(13,23),(13,21),(3,17),(3,19),(3,23),(3,21),(13,-17),(13,-19),(13,-23),(13,-21),(3,-17),(3,-19),(3,-23),(3,-21),(9,17),(9,19),(9,23),(9,21),(7,17),(7,19),(7,23),(7,21),(9,-17),(9,-19),(9,-23),(9,-21),(7,-17),(7,-19),(7,-23),(7,-21),(11,17),(11,19),(11,23),(11,21),(5,17),(5,19),(5,23),(5,21),(11,-17),(11,-19),(11,-23),(11,-21),(5,-17),(5,-19),(5,-23),(5,-21),(17,15),(17,13),(17,9),(17,11),(17,1),(17,3),(17,7),(17,5),(17,-15),(17,-13),(17,-9),(17,-11),(17,-1),(17,-3),(17,-7),(17,-5),(19,15),(19,13),(19,9),(19,11),(19,1),(19,3),(19,7),(19,5),(19,-15),(19,-13),(19,-9),(19,-11),(19,-1),(19,-3),(19,-7),(19,-5),(23,15),(23,13),(23,9),(23,11),(23,1),(23,3),(23,7),(23,5),(23,-15),(23,-13),(23,-9),(23,-11),(23,-1),(23,-3),(23,-7),(23,-5),(21,15),(21,13),(21,9),(21,11),(21,1),(21,3),(21,7),(21,5),(21,-15),(21,-13),(21,-9),(21,-11),(21,-1),(21,-3),(21,-7),(21,-5),(1,15),(1,13),(1,9),(1,11),(1,1),(1,3),(1,7),(1,5),(1,-15),(1,-13),(1,-9),(1,-11),(1,-1),(1,-3),(1,-7),(1,-5),(3,15),(3,13),(3,9),(3,11),(3,1),(3,3),(3,7),(3,5),(3,-15),(3,-13),(3,-9),(3,-11),(3,-1),(3,-3),(3,-7),(3,-5),(7,15),(7,13),(7,9),(7,11),(7,1),(7,3),(7,7),(7,5),(7,-15),(7,-13),(7,-9),(7,-11),(7,-1),(7,-3),(7,-7),(7,-5),(5,15),(5,13),(5,9),(5,11),(5,1),(5,3),(5,7),(5,5),(5,-15),(5,-13),(5,-9),(5,-11),(5,-1),(5,-3),(5,-7),(5,-5),(15,15),(15,13),(15,9),(15,11),(15,1),(15,3),(15,7),(15,5),(15,-15),(15,-13),(15,-9),(15,-11),(15,-1),(15,-3),(15,-7),(15,-5),(13,15),(13,13),(13,9),(13,11),(13,1),(13,3),(13,7),(13,5),(13,-15),(13,-13),(13,-9),(13,-11),(13,-1),(13,-3),(13,-7),(13,-5),(9,15),(9,13),(9,9),(9,11),(9,1),(9,3),(9,7),(9,5),(9,-15),(9,-13),(9,-9),(9,-11),(9,-1),(9,-3),(9,-7),(9,-5),(11,15),(11,13),(11,9),(11,11),(11,1),(11,3),(11,7),(11,5),(11,-15),(11,-13),(11,-9),(11,-11),(11,-1),(11,-3),(11,-7),(11,-5)]
const _QAM_CROSS_512_BIN = [(-15,17),(-15,19),(-15,21),(-15,23),(-1,23),(-1,21),(-1,19),(-1,17),(-1,-17),(-1,-19),(-1,-21),(-1,-23),(-15,-23),(-15,-21),(-15,-19),(-15,-17),(-13,17),(-13,19),(-13,21),(-13,23),(-3,23),(-3,21),(-3,19),(-3,17),(-3,-17),(-3,-19),(-3,-21),(-3,-23),(-13,-23),(-13,-21),(-13,-19),(-13,-17),(-11,17),(-11,19),(-11,21),(-11,23),(-5,23),(-5,21),(-5,19),(-5,17),(-5,-17),(-5,-19),(-5,-21),(-5,-23),(-11,-23),(-11,-21),(-11,-19),(-11,-17),(-9,17),(-9,19),(-9,21),(-9,23),(-7,23),(-7,21),(-7,19),(-7,17),(-7,-17),(-7,-19),(-7,-21),(-7,-23),(-9,-23),(-9,-21),(-9,-19),(-9,-17),(-23,15),(-23,13),(-23,11),(-23,9),(-23,7),(-23,5),(-23,3),(-23,1),(-23,-1),(-23,-3),(-23,-5),(-23,-7),(-23,-9),(-23,-11),(-23,-13),(-23,-15),(-21,15),(-21,13),(-21,11),(-21,9),(-21,7),(-21,5),(-21,3),(-21,1),(-21,-1),(-21,-3),(-21,-5),(-21,-7),(-21,-9),(-21,-11),(-21,-13),(-21,-15),(-19,15),(-19,13),(-19,11),(-19,9),(-19,7),(-19,5),(-19,3),(-19,1),(-19,-1),(-19,-3),(-19,-5),(-19,-7),(-19,-9),(-19,-11),(-19,-13),(-19,-15),(-17,15),(-17,13),(-17,11),(-17,9),(-17,7),(-17,5),(-17,3),(-17,1),(-17,-1),(-17,-3),(-17,-5),(-17,-7),(-17,-9),(-17,-11),(-17,-13),(-17,-15),(-15,15),(-15,13),(-15,11),(-15,9),(-15,7),(-15,5),(-15,3),(-15,1),(-15,-1),(-15,-3),(-15,-5),(-15,-7),(-15,-9),(-15,-11),(-15,-13),(-15,-15),(-13,15),(-13,13),(-13,11),(-13,9),(-13,7),(-13,5),(-13,3),(-13,1),(-13,-1),(-13,-3),(-13,-5),(-13,-7),(-13,-9),(-13,-11),(-13,-13),(-13,-15),(-11,15),(-11,13),(-11,11),(-11,9),(-11,7),(-11,5),(-11,3),(-11,1),(-11,-1),(-11,-3),(-11,-5),(-11,-7),(-11,-9),(-11,-11),(-11,-13),(-11,-15),(-9,15),(-9,13),(-9,11),(-9,9),(-9,7),(-9,5),(-9,3),(-9,1),(-9,-1),(-9,-3),(-9,-5),(-9,-7),(-9,-9),(-9,-11),(-9,-13),(-9,-15),(-7,15),(-7,13),(-7,11),(-7,9),(-7,7),(-7,5),(-7,3),(-7,1),(-7,-1),(-7,-3),(-7,-5),(-7,-7),(-7,-9),(-7,-11),(-7,-13),(-7,-15),(-5,15),(-5,13),(-5,11),(-5,9),(-5,7),(-5,5),(-5,3),(-5,1),(-5,-1),(-5,-3),(-5,-5),(-5,-7),(-5,-9),(-5,-11),(-5,-13),(-5,-15),(-3,15),(-3,13),(-3,11),(-3,9),(-3,7),(-3,5),(-3,3),(-3,1),(-3,-1),(-3,-3),(-3,-5),(-3,-7),(-3,-9),(-3,-11),(-3,-13),(-3,-15),(-1,15),(-1,13),(-1,11),(-1,9),(-1,7),(-1,5),(-1,3),(-1,1),(-1,-1),(-1,-3),(-1,-5),(-1,-7),(-1,-9),(-1,-11),(-1,-13),(-1,-15),(1,15),(1,13),(1,11),(1,9),(1,7),(1,5),(1,3),(1,1),(1,-1),(1,-3),(1,-5),(1,-7),(1,-9),(1,-11),(1,-13),(1,-15),(3,15),(3,13),(3,11),(3,9),(3,7),(3,5),(3,3),(3,1),(3,-1),(3,-3),(3,-5),(3,-7),(3,-9),(3,-11),(3,-13),(3,-15),(5,15),(5,13),(5,11),(5,9),(5,7),(5,5),(5,3),(5,1),(5,-1),(5,-3),(5,-5),(5,-7),(5,-9),(5,-11),(5,-13),(5,-15),(7,15),(7,13),(7,11),(7,9),(7,7),(7,5),(7,3),(7,1),(7,-1),(7,-3),(7,-5),(7,-7),(7,-9),(7,-11),(7,-13),(7,-15),(9,15),(9,13),(9,11),(9,9),(9,7),(9,5),(9,3),(9,1),(9,-1),(9,-3),(9,-5),(9,-7),(9,-9),(9,-11),(9,-13),(9,-15),(11,15),(11,13),(11,11),(11,9),(11,7),(11,5),(11,3),(11,1),(11,-1),(11,-3),(11,-5),(11,-7),(11,-9),(11,-11),(11,-13),(11,-15),(13,15),(13,13),(13,11),(13,9),(13,7),(13,5),(13,3),(13,1),(13,-1),(13,-3),(13,-5),(13,-7),(13,-9),(13,-11),(13,-13),(13,-15),(15,15),(15,13),(15,11),(15,9),(15,7),(15,5),(15,3),(15,1),(15,-1),(15,-3),(15,-5),(15,-7),(15,-9),(15,-11),(15,-13),(15,-15),(17,15),(17,13),(17,11),(17,9),(17,7),(17,5),(17,3),(17,1),(17,-1),(17,-3),(17,-5),(17,-7),(17,-9),(17,-11),(17,-13),(17,-15),(19,15),(19,13),(19,11),(19,9),(19,7),(19,5),(19,3),(19,1),(19,-1),(19,-3),(19,-5),(19,-7),(19,-9),(19,-11),(19,-13),(19,-15),(21,15),(21,13),(21,11),(21,9),(21,7),(21,5),(21,3),(21,1),(21,-1),(21,-3),(21,-5),(21,-7),(21,-9),(21,-11),(21,-13),(21,-15),(23,15),(23,13),(23,11),(23,9),(23,7),(23,5),(23,3),(23,1),(23,-1),(23,-3),(23,-5),(23,-7),(23,-9),(23,-11),(23,-13),(23,-15),(9,17),(9,19),(9,21),(9,23),(7,23),(7,21),(7,19),(7,17),(7,-17),(7,-19),(7,-21),(7,-23),(9,-23),(9,-21),(9,-19),(9,-17),(11,17),(11,19),(11,21),(11,23),(5,23),(5,21),(5,19),(5,17),(5,-17),(5,-19),(5,-21),(5,-23),(11,-23),(11,-21),(11,-19),(11,-17),(13,17),(13,19),(13,21),(13,23),(3,23),(3,21),(3,19),(3,17),(3,-17),(3,-19),(3,-21),(3,-23),(13,-23),(13,-21),(13,-19),(13,-17),(15,17),(15,19),(15,21),(15,23),(1,23),(1,21),(1,19),(1,17),(1,-17),(1,-19),(1,-21),(1,-23),(15,-23),(15,-21),(15,-19),(15,-17)]

const _QAM_RECTANGULAR_8_BIN_TO_GRAY = _qam_nonsquare_build_bin_to_gray(_QAM_RECTANGULAR_8_BIN, _QAM_RECTANGULAR_8_GRAY)
const _QAM_CROSS_32_BIN_TO_GRAY = _qam_nonsquare_build_bin_to_gray(_QAM_CROSS_32_BIN, _QAM_CROSS_32_GRAY)
const _QAM_CROSS_128_BIN_TO_GRAY = _qam_nonsquare_build_bin_to_gray(_QAM_CROSS_128_BIN, _QAM_CROSS_128_GRAY)
const _QAM_CROSS_512_BIN_TO_GRAY = _qam_nonsquare_build_bin_to_gray(_QAM_CROSS_512_BIN, _QAM_CROSS_512_GRAY)

const _QAM_RECTANGULAR_8_POSITION_INDEX = _qam_nonsquare_build_position_index(_QAM_RECTANGULAR_8_BIN, 4, 2)
const _QAM_CROSS_32_POSITION_INDEX = _qam_nonsquare_build_position_index(_QAM_CROSS_32_BIN, 6, 6)
const _QAM_CROSS_128_POSITION_INDEX = _qam_nonsquare_build_position_index(_QAM_CROSS_128_BIN, 12, 12)
const _QAM_CROSS_512_POSITION_INDEX = _qam_nonsquare_build_position_index(_QAM_CROSS_512_BIN, 24, 24)
