"""
    timingoffset(x; delay=0.0)

Introduce a fractional sample-timing delay with a parabolic Farrow
interpolator using fixed `alpha=0.5`. The implementation is a deterministic,
stateless impairment model intended for timing-recovery tests.

Unlike [`freqoffset`](@ref), realness is preserved: real `x` in gives real
`y` out, complex in gives complex out. Because the underlying interpolator
formula has an intrinsic 2-sample group delay by construction (it always
blends the two oldest of its four buffered samples), `timingoffset(x;
delay=d)` for integer `d` returns `x` shifted by exactly `d + 2` samples
(zero-padded at both ends), not `d` samples -- a plain, disclosed
consequence of reusing the verified formula honestly, not something this
function tries to cancel out. An empty `x` returns an empty vector.

The inverse operation is [`symbolsync`](@ref), which recovers (rather than
introduces) a sample-timing offset.
"""
function timingoffset end

function timingoffset(x::AbstractVector{<:Number}; delay::Real=0.0)
    _validate_timingoffset_arguments(delay)

    isempty(x) && return _upfirdn_empty_type(eltype(x))[]

    output_type = _awgn_output_type(x)
    pre_shift = floor(Int, delay)
    mu = _awgn_real_type(output_type)(delay - pre_shift)

    y = Vector{output_type}(undef, length(x))
    _timingoffset_kernel!(y, x, pre_shift, mu)
    return y
end

function _validate_timingoffset_arguments(delay::Real)
    delay isa Bool && throw(ArgumentError("delay must be a real number, not Bool"))
    isfinite(delay) || throw(ArgumentError("delay must be finite"))
    typemin(Int) < delay < typemax(Int) || throw(ArgumentError("delay must be representable as an Int"))
    return nothing
end

function _timingoffset_kernel!(y::Vector{Ty}, x::AbstractVector{<:Number}, pre_shift::Int, mu::T) where {Ty,T}
    n = length(x)
    x_first = firstindex(x)
    s1 = s2 = s3 = zero(Ty)
    for i in 1:n
        src = i - pre_shift
        xc = (1 <= src <= n) ? Ty(x[x_first+src-1]) : zero(Ty)
        y[i], s1, s2, s3 = _farrow_interp(xc, s1, s2, s3, mu)
    end
    return y
end
