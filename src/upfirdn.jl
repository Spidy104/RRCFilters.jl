"""
    upfirdn(x, h, p=1, q=1)

Upsample `x` by `p`, filter with the FIR coefficients `h`, and downsample by
`q`: zero-stuffed upsampling with no
trailing padding after the last sample, full convolution with `h`, then
decimation starting at the first sample. `p` and `q` default to `1` and must
be positive integers, not `Bool`.

The result length is `ceil(((length(x) - 1) * p + length(h)) / q)`. An empty
`x` returns an empty vector; `h` must be nonempty, and both `x` and `h` must
contain only finite values. Concrete floating-point element types are
preserved and precision is promoted across `x` and `h` the same way `awgn`
promotes signal precision; integer and rational elements promote to
`Float64`. The result is a newly allocated vector, and neither `x` nor `h`
is modified.
"""
function upfirdn end

function upfirdn(x::AbstractVector{<:Number}, h::AbstractVector{<:Number}, p::Integer=1, q::Integer=1)
    nx = length(x)
    nh = length(h)
    _validate_upfirdn_arguments(nx, nh, p, q)
    output_type, h_type = _upfirdn_output_type(x, h)
    isempty(x) && return output_type[]
    nx == 1 && return _upfirdn_single_sample(x, h, Int(q), nh, output_type, h_type)

    return _upfirdn_kernel(x, h, Int(p), Int(q), nx, nh, output_type, h_type)
end

function _upfirdn_single_sample(x::AbstractVector{<:Number}, h::AbstractVector{<:Number}, q::Int, nh::Int,
                                ::Type{T}, ::Type{H}) where {T,H}
    y = Vector{T}(undef, cld(nh, q))
    sample = T(x[firstindex(x)])
    h_first = firstindex(h)
    for output_index in eachindex(y)
        tap_offset = (output_index - 1) * q
        y[output_index] = sample * T(H(h[h_first + tap_offset]))
    end
    return y
end

function _upfirdn_phase_shape(p_int::Int)
    is_pow2 = ispow2(p_int)
    shift = is_pow2 ? trailing_zeros(p_int) : 0
    mask = p_int - 1
    return is_pow2, shift, mask
end

@inline function _upfirdn_quotient_phase(mm1::Int, p_int::Int, is_pow2::Bool, shift::Int, mask::Int)
    p_int == 1 && return mm1, 0
    is_pow2 && return mm1 >> shift, mm1 & mask
    return fldmod(mm1, p_int)
end

# Output-index range [lo, hi] over which EVERY phase's tap window is fully
# inside x: no start-of-signal truncation (unclamped_top >= phase_length)
# and no end-of-signal truncation (unclamped_top <= nx). Inside this range
# the generic loop's skip/i_end bounds (see below) always take their
# trivial values (0 and phase_length), so a dedicated branch-free loop can
# replace the per-sample min() bookkeeping entirely. Derived in closed form
# only for the two shapes this library's own upfirdn calls ever use --
# pure upsampling (q_int==1, any p_int) and pure decimation (p_int==1, any
# q_int); general rational resampling (both >1) falls back to "no fast
# region" (lo > hi), which keeps every output on the existing safe
# per-sample loop, unchanged from before this optimization.
function _upfirdn_steady_region(nx::Int, max_phase_length::Int, p_int::Int, q_int::Int, output_length::Int)
    lo_quotient = max_phase_length - 1
    hi_quotient = nx - 1
    lo_quotient > hi_quotient && return output_length + 1, output_length
    q_int == 1 && return lo_quotient * p_int + 1, min(hi_quotient * p_int + p_int, output_length)
    p_int == 1 && return cld(lo_quotient, q_int) + 1, min(fld(hi_quotient, q_int) + 1, output_length)
    return output_length + 1, output_length
end

function _upfirdn_boundary_range!(y::Vector{T}, x::AbstractVector{<:Number}, phase_data::Vector{H}, phase_offsets::Vector{Int},
                                   x_first::Int, p_int::Int, q_int::Int, is_pow2::Bool, shift::Int, mask::Int,
                                   nx::Int, from::Int, to::Int) where {T,H}
    @inbounds for output_index in from:to
        m = (output_index - 1) * q_int + 1
        quotient, phase = _upfirdn_quotient_phase(m - 1, p_int, is_pow2, shift, mask)
        unclamped_top = quotient + 1
        n_top = min(nx, unclamped_top)
        skip = unclamped_top - n_top
        phase_start = phase_offsets[phase + 1]
        phase_length = phase_offsets[phase + 2] - phase_start
        i_end = min(phase_length, unclamped_top)
        accumulator = zero(T)
        @simd for i in (skip + 1):i_end
            accumulator += phase_data[phase_start + i] * T(x[x_first + unclamped_top - i])
        end
        y[output_index] = accumulator
    end
    return nothing
end

function _upfirdn_kernel(x::AbstractVector{<:Number}, h::AbstractVector{<:Number}, p_int::Int, q_int::Int, nx::Int, nh::Int, ::Type{T}, ::Type{H}) where {T,H}
    output_length = _upfirdn_output_length(nx, nh, p_int, q_int)
    y = Vector{T}(undef, output_length)

    x_first = firstindex(x)
    phase_data, phase_offsets = _upfirdn_polyphase_filters(h, firstindex(h), nh, p_int, H)
    is_pow2, shift, mask = _upfirdn_phase_shape(p_int)
    max_phase_length = maximum(phase_offsets[phase + 2] - phase_offsets[phase + 1] for phase in 0:(p_int - 1))
    lo_out, hi_out = _upfirdn_steady_region(nx, max_phase_length, p_int, q_int, output_length)

    _upfirdn_boundary_range!(y, x, phase_data, phase_offsets, x_first, p_int, q_int, is_pow2, shift, mask, nx, 1, lo_out - 1)

    @inbounds for output_index in lo_out:hi_out
        m = (output_index - 1) * q_int + 1
        quotient, phase = _upfirdn_quotient_phase(m - 1, p_int, is_pow2, shift, mask)
        unclamped_top = quotient + 1
        phase_start = phase_offsets[phase + 1]
        phase_length = phase_offsets[phase + 2] - phase_start
        accumulator = zero(T)
        @simd for i in 1:phase_length
            accumulator += phase_data[phase_start + i] * T(x[x_first + unclamped_top - i])
        end
        y[output_index] = accumulator
    end

    _upfirdn_boundary_range!(y, x, phase_data, phase_offsets, x_first, p_int, q_int, is_pow2, shift, mask, nx, hi_out + 1, output_length)

    return y
end

const _UpfirdnSIMDFloat = Union{Float32,Float64}

function _upfirdn_complex_boundary_range!(y::Vector{Complex{R}}, x::AbstractVector{<:Number}, phase_data::Vector{H}, phase_offsets::Vector{Int},
                                           x_first::Int, p_int::Int, q_int::Int, is_pow2::Bool, shift::Int, mask::Int,
                                           nx::Int, from::Int, to::Int) where {R,H}
    @inbounds for output_index in from:to
        m = (output_index - 1) * q_int + 1
        quotient, phase = _upfirdn_quotient_phase(m - 1, p_int, is_pow2, shift, mask)
        unclamped_top = quotient + 1
        n_top = min(nx, unclamped_top)
        skip = unclamped_top - n_top
        phase_start = phase_offsets[phase + 1]
        phase_length = phase_offsets[phase + 2] - phase_start
        i_end = min(phase_length, unclamped_top)
        acc_re = zero(R)
        acc_im = zero(R)
        @simd for i in (skip + 1):i_end
            tap = phase_data[phase_start + i]
            sample = x[x_first + unclamped_top - i]
            acc_re += tap * real(sample)
            acc_im += tap * imag(sample)
        end
        y[output_index] = Complex{R}(acc_re, acc_im)
    end
    return nothing
end

# Julia's generic Complex(::Real, ::Real)-accumulator loop
# (`accumulator += tap * Complex{R}(x[...])`) does not auto-vectorize --
# confirmed via @code_native, no ymm/zmm/fma emitted regardless of tap
# count -- because the compiler doesn't treat a packed Complex{R} reduction
# as a vectorizable pattern the way it does two independent real
# accumulators. Splitting the accumulator into separate real/imaginary
# scalars (still reading directly from the interleaved x, no upfront
# array copy needed) restores full-width AVX2/AVX-512+FMA vectorization of
# the reduction, and, measured directly at 1,000,000-symbol scale, never
# loses to the un-vectorized combined-accumulator form -- a wash for very
# short per-phase tap counts, a large win (1.6x-4.6x here) for medium/long
# ones -- so unlike the previous (now-replaced) full-array-copy SoA
# approach, this needs no tap-count threshold: reading directly from x
# instead of a copied buffer also removes that approach's O(nx) upfront
# copy cost, which used to dominate for heavily-decimated receive
# filtering (nx greatly exceeding the output length).
const _UPFIRDN_INTERP_MAX_P = 16

function _upfirdn_interp_taps(h::AbstractVector{<:Number}, h_first::Int, nh::Int, p::Int, ::Type{H}) where {H}
    L = cld(nh, p)
    t = zeros(H, p, L)
    for phase in 0:(p - 1), i in 1:L
        idx = phase + (i - 1) * p
        idx < nh && (t[phase + 1, i] = H(h[h_first + idx]))
    end
    return t, L
end

@generated function _upfirdn_interp_steady!(y::Vector{Complex{R}}, x::AbstractVector{<:Number}, t::Matrix{H},
                                             L::Int, lo_q::Int, hi_q::Int, x_first::Int, ::Val{P}) where {R,H,P}
    quote
        @inbounds for q in lo_q:hi_q
            top = q + 1
            Base.Cartesian.@nexprs $P k -> (ar_k = zero(R); ai_k = zero(R))
            for i in 1:L
                s = x[x_first + top - i]
                re = R(real(s))
                im = R(imag(s))
                Base.Cartesian.@nexprs $P k -> begin
                    c_k = t[k, i]
                    ar_k = muladd(c_k, re, ar_k)
                    ai_k = muladd(c_k, im, ai_k)
                end
            end
            base = q * $P
            Base.Cartesian.@nexprs $P k -> (y[base + k] = Complex{R}(ar_k, ai_k))
        end
        return nothing
    end
end

const _UPFIRDN_DECIM_UNROLL_MAX_TAPS = 64

function _upfirdn_decim_steady!(y::Vector{Complex{R}}, x::AbstractVector{<:Number}, taps::Vector{H},
                                 nh::Int, q_int::Int, x_first::Int, lo_m::Int, hi_m::Int) where {R,H}
    if nh > _UPFIRDN_DECIM_UNROLL_MAX_TAPS
        @inbounds for m in lo_m:hi_m
            top = (m - 1) * q_int + 1
            ar = zero(R)
            ai = zero(R)
            @simd for i in 1:nh
                s = x[x_first + top - i]
                ar += taps[i] * R(real(s))
                ai += taps[i] * R(imag(s))
            end
            y[m] = Complex{R}(ar, ai)
        end
        return nothing
    end
    m = lo_m
    @inbounds while m + 3 <= hi_m
        t1 = (m - 1) * q_int + 1
        t2 = t1 + q_int
        t3 = t2 + q_int
        t4 = t3 + q_int
        a1r = zero(R); a1i = zero(R); a2r = zero(R); a2i = zero(R)
        a3r = zero(R); a3i = zero(R); a4r = zero(R); a4i = zero(R)
        for i in 1:nh
            c = taps[i]
            s1 = x[x_first + t1 - i]
            s2 = x[x_first + t2 - i]
            s3 = x[x_first + t3 - i]
            s4 = x[x_first + t4 - i]
            a1r = muladd(c, R(real(s1)), a1r); a1i = muladd(c, R(imag(s1)), a1i)
            a2r = muladd(c, R(real(s2)), a2r); a2i = muladd(c, R(imag(s2)), a2i)
            a3r = muladd(c, R(real(s3)), a3r); a3i = muladd(c, R(imag(s3)), a3i)
            a4r = muladd(c, R(real(s4)), a4r); a4i = muladd(c, R(imag(s4)), a4i)
        end
        y[m] = Complex{R}(a1r, a1i)
        y[m + 1] = Complex{R}(a2r, a2i)
        y[m + 2] = Complex{R}(a3r, a3i)
        y[m + 3] = Complex{R}(a4r, a4i)
        m += 4
    end
    @inbounds while m <= hi_m
        top = (m - 1) * q_int + 1
        ar = zero(R)
        ai = zero(R)
        @simd for i in 1:nh
            s = x[x_first + top - i]
            ar += taps[i] * R(real(s))
            ai += taps[i] * R(imag(s))
        end
        y[m] = Complex{R}(ar, ai)
        m += 1
    end
    return nothing
end

function _upfirdn_kernel(x::AbstractVector{<:Number}, h::AbstractVector{<:Number}, p_int::Int, q_int::Int, nx::Int, nh::Int, ::Type{Complex{R}}, ::Type{H}) where {R<:_UpfirdnSIMDFloat,H<:_UpfirdnSIMDFloat}
    output_length = _upfirdn_output_length(nx, nh, p_int, q_int)
    y = Vector{Complex{R}}(undef, output_length)
    phase_data, phase_offsets = _upfirdn_polyphase_filters(h, firstindex(h), nh, p_int, H)
    is_pow2, shift, mask = _upfirdn_phase_shape(p_int)
    x_first = firstindex(x)
    max_phase_length = maximum(phase_offsets[phase + 2] - phase_offsets[phase + 1] for phase in 0:(p_int - 1))
    lo_out, hi_out = _upfirdn_steady_region(nx, max_phase_length, p_int, q_int, output_length)

    _upfirdn_complex_boundary_range!(y, x, phase_data, phase_offsets, x_first, p_int, q_int, is_pow2, shift, mask, nx, 1, lo_out - 1)

    if q_int == 1 && 2 <= p_int <= _UPFIRDN_INTERP_MAX_P && nh >= p_int && lo_out <= hi_out
        t, L = _upfirdn_interp_taps(h, firstindex(h), nh, p_int, H)
        _upfirdn_interp_steady!(y, x, t, L, max_phase_length - 1, nx - 1, x_first, Val(p_int))
        _upfirdn_complex_boundary_range!(y, x, phase_data, phase_offsets, x_first, p_int, q_int, is_pow2, shift, mask, nx, hi_out + 1, output_length)
        return y
    end
    if p_int == 1
        _upfirdn_decim_steady!(y, x, phase_data, nh, q_int, x_first, lo_out, hi_out)
        _upfirdn_complex_boundary_range!(y, x, phase_data, phase_offsets, x_first, p_int, q_int, is_pow2, shift, mask, nx, hi_out + 1, output_length)
        return y
    end

    @inbounds for output_index in lo_out:hi_out
        m = (output_index - 1) * q_int + 1
        quotient, phase = _upfirdn_quotient_phase(m - 1, p_int, is_pow2, shift, mask)
        unclamped_top = quotient + 1
        phase_start = phase_offsets[phase + 1]
        phase_length = phase_offsets[phase + 2] - phase_start
        acc_re = zero(R)
        acc_im = zero(R)
        @simd for i in 1:phase_length
            tap = phase_data[phase_start + i]
            sample = x[x_first + unclamped_top - i]
            acc_re += tap * real(sample)
            acc_im += tap * imag(sample)
        end
        y[output_index] = Complex{R}(acc_re, acc_im)
    end

    _upfirdn_complex_boundary_range!(y, x, phase_data, phase_offsets, x_first, p_int, q_int, is_pow2, shift, mask, nx, hi_out + 1, output_length)

    return y
end

function _upfirdn_polyphase_filters(h::AbstractVector{<:Number}, h_first::Int, nh::Int, p::Int, ::Type{T}) where {T}
    offsets = Vector{Int}(undef, p + 1)
    offsets[1] = 0
    for phase in 0:(p - 1)
        offsets[phase + 2] = offsets[phase + 1] + max(0, fld(nh - phase - 1, p) + 1)
    end

    data = Vector{T}(undef, nh)
    for phase in 0:(p - 1)
        base = offsets[phase + 1]
        count = offsets[phase + 2] - base
        for i in 1:count
            data[base + i] = T(h[h_first + phase + (i - 1) * p])
        end
    end

    return data, offsets
end

function _validate_upfirdn_arguments(nx::Int, nh::Int, p::Integer, q::Integer)
    p isa Bool && throw(ArgumentError("p must be a positive integer, not Bool"))
    q isa Bool && throw(ArgumentError("q must be a positive integer, not Bool"))
    p > 0 || throw(ArgumentError("p must be positive"))
    q > 0 || throw(ArgumentError("q must be positive"))
    p <= typemax(Int) || throw(ArgumentError("p is too large"))
    q <= typemax(Int) || throw(ArgumentError("q is too large"))
    nh > 0 || throw(ArgumentError("h must be nonempty"))
    nx == 0 && return nothing
    (nx - 1) <= div(typemax(Int) - nh, Int(p)) || throw(ArgumentError("x, h, and p combination is too large"))
    return nothing
end

function _upfirdn_output_type(x::AbstractVector{<:Number}, h::AbstractVector{<:Number})
    x_type = isempty(x) ? _upfirdn_empty_type(eltype(x)) : _awgn_output_type(x, "x")
    h_type = _awgn_output_type(h, "h")
    real_type = promote_type(_awgn_real_type(x_type), _awgn_real_type(h_type))
    output_type = (x_type <: Complex || h_type <: Complex) ? Complex{real_type} : real_type
    return output_type, h_type
end

function _upfirdn_output_length(nx::Int, nh::Int, p::Int, q::Int)
    iszero(nx) && return 0
    return cld((nx - 1) * p + nh, q)
end
