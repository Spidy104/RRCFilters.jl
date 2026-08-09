"""
    rcosdesign(beta, span, sps; shape=:sqrt)

Design a unit-energy raised-cosine pulse-shaping FIR filter.

`beta` is the roll-off factor in `[0, 1]`; `span` is the filter length in
symbols; and `sps` is the number of samples per symbol. `span * sps` must be
even so the returned filter has an odd number of samples and a sample at zero.

Set `shape=:sqrt` for a root-raised-cosine filter or `shape=:normal` for a
raised-cosine filter.
"""
function rcosdesign(beta::Real, span::Integer, sps::Integer; shape::Symbol=:sqrt)
    beta isa Bool && throw(ArgumentError("beta must be a real number, not Bool"))
    span isa Bool && throw(ArgumentError("span must be an integer, not Bool"))
    sps isa Bool && throw(ArgumentError("sps must be an integer, not Bool"))
    isfinite(beta) || throw(ArgumentError("beta must be finite"))
    zero(beta) <= beta <= one(beta) || throw(ArgumentError("beta must be in [0, 1]"))
    span > 0 || throw(ArgumentError("span must be positive"))
    sps > 0 || throw(ArgumentError("sps must be positive"))
    span <= div(typemax(Int) - 1, sps) || throw(ArgumentError("span * sps is too large"))
    iseven(span * sps) || throw(ArgumentError("span * sps must be even"))
    shape in (:sqrt, :normal) || throw(ArgumentError("shape must be :sqrt or :normal"))
    T = beta isa AbstractFloat ? typeof(beta) : Float64
    b = T(beta)
    piT = T(pi)
    sample_count = Int(span * sps)
    half_count = sample_count ÷ 2
    coefficients = Vector{T}(undef, sample_count + 1)
    for (index, n) in enumerate(-half_count:half_count)
        t = T(n) / T(sps)
        coefficients[index] = shape === :sqrt ? _root_raised_cosine(t, b, piT) : _raised_cosine(t, b, piT)
    end
    coefficients ./= sqrt(sum(abs2, coefficients))
    return coefficients
end

function _raised_cosine(t::T, beta::T, piT::T) where {T<:AbstractFloat}
    iszero(beta) && return sinc(t)
    iszero(t) && return one(T)
    scaled_time = T(2) * beta * t
    if isapprox(abs(scaled_time), one(T); rtol=sqrt(eps(T)) / T(4), atol=sqrt(eps(T)) / T(4))
        return beta / T(2) * sin(piT / (T(2) * beta))
    end
    return sinc(t) * cos(piT * beta * t) / (one(T) - scaled_time^2)
end

function _root_raised_cosine(t::T, beta::T, piT::T) where {T<:AbstractFloat}
    iszero(beta) && return sinc(t)
    iszero(t) && return one(T) - beta + T(4) * beta / piT
    scaled_time = T(4) * beta * t
    if isapprox(abs(scaled_time), one(T); rtol=sqrt(eps(T)) / T(4), atol=sqrt(eps(T)) / T(4))
        angle = piT / (T(4) * beta)
        return beta / sqrt(T(2)) * ((one(T) + T(2) / piT) * sin(angle) + (one(T) - T(2) / piT) * cos(angle))
    end
    numerator = sin(piT * t * (one(T) - beta)) + T(4) * beta * t * cos(piT * t * (one(T) + beta))
    denominator = piT * t * (one(T) - scaled_time^2)
    return numerator / denominator
end
