_upfirdn_empty_type(::Type{T}) where {T<:AbstractFloat} = T
_upfirdn_empty_type(::Type{Complex{T}}) where {T<:AbstractFloat} = Complex{T}
_upfirdn_empty_type(::Type{<:Complex}) = ComplexF64
_upfirdn_empty_type(::Type{<:Real}) = Float64
_upfirdn_empty_type(::Type) = Float64

_awgn_real_type(::Type{T}) where {T<:AbstractFloat} = T
_awgn_real_type(::Type{Complex{T}}) where {T<:AbstractFloat} = T

# Fast path for the common case: x is already a concretely and homogeneously
# typed vector of one of the four floating-point types this library ever
# promotes to (or their Complex wrapping), so there is nothing to promote
# and every element is already known to be Real/Complex -- the only
# per-element work actually needed is the NaN/Inf check, which `all` can
# loop over with a plain vectorizable predicate instead of the generic
# method's per-sample type-dispatch/promotion bookkeeping (which showed up
# as a real cost -- several thousand microseconds at 1,000,000+ elements,
# dominated by GC pauses from the boxed `Union{Nothing,Type}` running
# `float_type` -- while investigating upfirdn's wrapper overhead). The type
# union below is deliberately an explicit enumeration of concrete leaf
# types, NOT the abstract `AbstractFloat`/`Complex{<:AbstractFloat}`
# supertypes -- an abstractly-typed vector (e.g. `AbstractFloat[...]`, which
# can mix Float32 and BigFloat elements at runtime) must still take the slow
# path so mixed-precision promotion is computed correctly; dispatching on
# the concrete eltype `T` of the array is what guarantees homogeneity here.
const _AwgnConcreteFloatType = Union{
    Float16,Float32,Float64,BigFloat,
    Complex{Float16},Complex{Float32},Complex{Float64},Complex{BigFloat},
}

function _awgn_output_type(x::AbstractVector{T}, name::AbstractString="x") where {T<:_AwgnConcreteFloatType}
    all(isfinite, x) || throw(ArgumentError("$name must contain only finite samples"))
    return T
end

function _awgn_output_type(x::AbstractVector{<:Number}, name::AbstractString="x")
    float_type = nothing
    complex_signal = false

    for sample in x
        sample isa Union{Real,Complex} || throw(ArgumentError("$name must contain real or complex samples"))
        isfinite(sample) || throw(ArgumentError("$name must contain only finite samples"))
        component = sample isa Complex ? real(sample) : sample
        component isa Real || throw(ArgumentError("$name must contain real or complex samples"))
        complex_signal |= sample isa Complex

        if component isa AbstractFloat
            float_type = isnothing(float_type) ? typeof(component) : promote_type(float_type, typeof(component))
        elseif !isnothing(float_type)
            float_type = promote_type(float_type, Float64)
        end
    end

    real_type = isnothing(float_type) ? Float64 : float_type
    real_type <: Union{Float16,Float32,Float64,BigFloat} || throw(ArgumentError("$name uses an unsupported floating-point type"))
    return complex_signal ? Complex{real_type} : real_type
end

function _farrow_interp(xc::Ty, s1::Ty, s2::Ty, s3::Ty, mu::T) where {Ty<:Number,T<:AbstractFloat}
    c0 = s2
    c1 = -T(0.5) * xc + T(1.5) * s1 - T(0.5) * s2 - T(0.5) * s3
    c2 = T(0.5) * xc - T(0.5) * s1 - T(0.5) * s2 + T(0.5) * s3
    y = c0 + mu * c1 + mu * mu * c2
    return y, xc, s1, s2
end

_binary_to_gray(value::Integer) = value ⊻ (value >> 1)

function _gray_to_binary(value::Integer)
    binary = zero(value)
    while !iszero(value)
        binary ⊻= value
        value >>= 1
    end
    return binary
end

function _validate_bit_vector(bits::AbstractVector{<:Integer})
    for bit in bits
        (bit == 0 || bit == 1) || throw(ArgumentError("bits must contain only 0 or 1"))
    end
end
function _bits_to_integer(bits::AbstractVector{<:Integer}, start::Int, width::Int)
    value = 0
    for index in start:(start + width - 1)
        value = (value << 1) | Int(bits[index])
    end
    return value
end

function _integer_to_bits!(bits::BitVector, start::Int, width::Int, value::Int)
    for offset in 0:(width - 1)
        bits[start + offset] = !iszero(value & (1 << (width - offset - 1)))
    end
end
