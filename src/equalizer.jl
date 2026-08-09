const _EQUALIZE_DEFAULT_CONSTELLATION = ComplexF64[(1 + im), (-1 + im), (-1 - im), (1 - im)] ./ sqrt(2.0)

"""
    equalize(x, training_symbols=ComplexF64[]; algorithm=:lms, num_forward_taps=5,
             num_feedback_taps=0, reference_tap=3, constellation=<unit-power QPSK>,
             step_size=0.01, forgetting_factor=0.99, initial_inverse_correlation=0.1,
             projection_order=4, regularization=1.0)

Adaptive linear or decision-feedback equalizer. Returns
`(equalized, error, weights)`:
`equalized`/`error` have length `length(x)`; `weights` always has length
`num_forward_taps + num_feedback_taps` regardless of `length(x)` (including
when `x` is empty), since it is configuration-shaped, not data-length-shaped.
Element type is `Complex{T}` if any of `x`/`constellation`/(nonempty)
`training_symbols` is complex, otherwise plain `T`; realness is preserved.

`algorithm` is `:lms` (fixed step size), `:rls` (recursive least squares,
forgetting factor), or `:apa` (affine projection). Decision feedback extends
the combined tap-delay line rather than selecting a separate algorithm.
The combined
tap-delay line `u` has length `num_forward_taps + num_feedback_taps`: the
first `num_forward_taps` slots hold input samples (`u[1]` = newest, shifted
in every step, same convention as before), the remaining `num_feedback_taps`
slots hold past `ref` values (`u[num_forward_taps + 1]` = most recent `ref`,
one step delayed -- i.e. the feedback line is fed by whatever `ref` was
computed on the *previous* step, not a freshly recomputed decision). Output
`y = sum_k conj(w[k]) * u[k]` is one dot product over the whole combined
vector; error `e = ref - y`. `ref` is the training symbol delay-compensated
by `reference_tap - 1` samples while training symbols remain, a hard
decision (nearest `constellation` point, ties broken toward the
first/lowest-index point) during the initial `reference_tap - 1`-sample
warm-up (no valid delayed training symbol exists yet) and after training
symbols are exhausted (decision-directed operation) -- unchanged by
`num_feedback_taps`. `:lms` updates `w[k] += step_size * u[k] * conj(e)` for
every tap in the combined vector; `:rls` maintains an inverse-correlation
matrix `P` (init `initial_inverse_correlation * I`, sized
`(num_forward_taps + num_feedback_taps)^2`) and forgetting factor
`forgetting_factor`. Weights start at zero.

`num_feedback_taps=0` selects linear-only operation; positive values enable
decision feedback. `num_forward_taps` must be at least one.

`:apa` implements the standard affine projection algorithm. Per step it
maintains a sliding window of the last `projection_order` combined tap
vectors `u` (`Uhist`, shifted like `u` itself) and reference values (`dhist`);
solves the small `projection_order`-square linear system
`R * a = conj.(dhist .- (weights'Uhist))` (`R` the window's Gram matrix plus
`regularization` on the diagonal, recomputed fresh every step rather than
incrementally maintained) via in-place Gauss-Jordan elimination with partial
pivoting, then
`weights .+= step_size .* (Uhist * a)`. `regularization` (default `1.0`)
stabilizes the Gram matrix. `step_size` must be in `(0, 1]` for `:apa` (normalized step
sizes are only meaningful in that range) -- `:lms`'s own `step_size` keeps
its unbounded-above `> 0` requirement unchanged. `projection_order >= 1` is
permitted; `projection_order=1` is exactly NLMS.

Empty `training_symbols` enables pure decision-directed operation from
zero-initialized weights. `length(training_symbols) <= length(x)` and
`length(constellation) >= 2` are required. Out of scope: `:cma` blind
equalization, fractional/oversampled (`InputSamplesPerSymbol`) operation,
`InputDelay`, and non-scalar `InitialInverseCorrelationMatrix`.

An empty `x` returns empty `equalized`/`error` and an all-zero `weights` of
length `num_forward_taps + num_feedback_taps`.
"""
function equalize end

function equalize(x::AbstractVector{<:Number}, training_symbols::AbstractVector{<:Number}=ComplexF64[];
                   algorithm::Symbol=:lms, num_forward_taps::Integer=5, num_feedback_taps::Integer=0,
                   reference_tap::Integer=3,
                   constellation::AbstractVector{<:Number}=_EQUALIZE_DEFAULT_CONSTELLATION,
                   step_size::Real=0.01, forgetting_factor::Real=0.99,
                   initial_inverse_correlation::Real=0.1,
                   projection_order::Integer=4, regularization::Real=1.0)
    _validate_equalize_arguments(length(training_symbols), length(x), algorithm, num_forward_taps,
                                  num_feedback_taps, reference_tap, step_size, forgetting_factor,
                                  initial_inverse_correlation, projection_order, regularization,
                                  constellation)
    num_taps = Int(num_forward_taps) + Int(num_feedback_taps)

    if isempty(x)
        x_empty_type = _upfirdn_empty_type(eltype(x))
        constellation_type = _awgn_output_type(constellation, "constellation")
        real_type = promote_type(_awgn_real_type(x_empty_type), _awgn_real_type(constellation_type))
        output_type = (x_empty_type <: Complex || constellation_type <: Complex) ? Complex{real_type} : real_type
        return output_type[], output_type[], zeros(output_type, num_taps)
    end

    output_type = _equalize_output_type(x, constellation, training_symbols)
    equalized = Vector{output_type}(undef, length(x))
    err = Vector{output_type}(undef, length(x))
    weights = zeros(output_type, num_taps)

    if algorithm === :lms
        _equalize_kernel!(equalized, err, weights, x, training_symbols, constellation,
                           Int(reference_tap), Int(num_forward_taps), step_size, Val(:lms))
    elseif algorithm === :rls
        _equalize_kernel!(equalized, err, weights, x, training_symbols, constellation,
                           Int(reference_tap), Int(num_forward_taps), forgetting_factor,
                           initial_inverse_correlation, Val(:rls))
    else
        _equalize_kernel!(equalized, err, weights, x, training_symbols, constellation,
                           Int(reference_tap), Int(num_forward_taps), Int(projection_order),
                           step_size, regularization, Val(:apa))
    end
    return equalized, err, weights
end

function _validate_equalize_arguments(training_length::Int, x_length::Int, algorithm::Symbol,
                                       num_forward_taps::Integer, num_feedback_taps::Integer,
                                       reference_tap::Integer, step_size::Real, forgetting_factor::Real,
                                       initial_inverse_correlation::Real, projection_order::Integer,
                                       regularization::Real, constellation::AbstractVector{<:Number})
    training_length <= x_length || throw(ArgumentError("training_symbols must not be longer than x"))

    algorithm in (:lms, :rls, :apa) || throw(ArgumentError("algorithm must be :lms, :rls, or :apa"))

    num_forward_taps isa Bool && throw(ArgumentError("num_forward_taps must be an integer, not Bool"))
    num_forward_taps >= 1 || throw(ArgumentError("num_forward_taps must be >= 1"))

    num_feedback_taps isa Bool && throw(ArgumentError("num_feedback_taps must be an integer, not Bool"))
    num_feedback_taps >= 0 || throw(ArgumentError("num_feedback_taps must be >= 0"))

    reference_tap isa Bool && throw(ArgumentError("reference_tap must be an integer, not Bool"))
    (reference_tap >= 1 && reference_tap <= num_forward_taps) ||
        throw(ArgumentError("reference_tap must be in [1, num_forward_taps]"))

    projection_order isa Bool && throw(ArgumentError("projection_order must be an integer, not Bool"))
    projection_order >= 1 || throw(ArgumentError("projection_order must be >= 1"))

    step_size isa Bool && throw(ArgumentError("step_size must be a real number, not Bool"))
    isfinite(step_size) || throw(ArgumentError("step_size must be finite"))
    step_size > 0 || throw(ArgumentError("step_size must be positive"))
    (algorithm === :apa && step_size > 1) &&
        throw(ArgumentError("step_size must be in (0, 1] for algorithm=:apa"))

    forgetting_factor isa Bool && throw(ArgumentError("forgetting_factor must be a real number, not Bool"))
    isfinite(forgetting_factor) || throw(ArgumentError("forgetting_factor must be finite"))
    (forgetting_factor > 0 && forgetting_factor <= 1) || throw(ArgumentError("forgetting_factor must be in (0, 1]"))

    initial_inverse_correlation isa Bool && throw(ArgumentError("initial_inverse_correlation must be a real number, not Bool"))
    isfinite(initial_inverse_correlation) || throw(ArgumentError("initial_inverse_correlation must be finite"))
    initial_inverse_correlation > 0 || throw(ArgumentError("initial_inverse_correlation must be positive"))

    regularization isa Bool && throw(ArgumentError("regularization must be a real number, not Bool"))
    isfinite(regularization) || throw(ArgumentError("regularization must be finite"))
    regularization > 0 || throw(ArgumentError("regularization must be positive"))

    length(constellation) >= 2 || throw(ArgumentError("constellation must have at least 2 points"))

    return nothing
end

function _equalize_output_type(x::AbstractVector{<:Number}, constellation::AbstractVector{<:Number},
                                training_symbols::AbstractVector{<:Number})
    t_x = _awgn_output_type(x, "x")
    t_c = _awgn_output_type(constellation, "constellation")
    t_train = isempty(training_symbols) ? t_c : _awgn_output_type(training_symbols, "training_symbols")
    real_type = promote_type(_awgn_real_type(t_x), _awgn_real_type(t_c), _awgn_real_type(t_train))
    is_complex = (t_x <: Complex) || (t_c <: Complex) || (t_train <: Complex)
    return is_complex ? Complex{real_type} : real_type
end

function _equalizer_shift!(u::AbstractVector{T}, newest::T) where {T}
    n = length(u)
    @inbounds for k in n:-1:2
        u[k] = u[k - 1]
    end
    @inbounds u[1] = newest
    return u
end

function _equalizer_output(weights::Vector{T}, u::Vector{T}) where {T}
    y = zero(T)
    for k in eachindex(weights, u)
        y += conj(weights[k]) * u[k]
    end
    return y
end

function _equalizer_shift_matrix_columns!(U::Matrix{T}, newest::AbstractVector{T}) where {T}
    L, N = size(U)
    @inbounds for j in N:-1:2, i in 1:L
        U[i, j] = U[i, j - 1]
    end
    @inbounds for i in 1:L
        U[i, 1] = newest[i]
    end
    return U
end

function _equalizer_apa_solve!(R::Matrix{T}, evec::Vector{T}) where {T}
    n = length(evec)
    for k in 1:n
        pivot_row = k
        pivot_mag = abs(R[k, k])
        for i in (k + 1):n
            mag = abs(R[i, k])
            if mag > pivot_mag
                pivot_mag = mag
                pivot_row = i
            end
        end
        if pivot_row != k
            for j in 1:n
                R[k, j], R[pivot_row, j] = R[pivot_row, j], R[k, j]
            end
            evec[k], evec[pivot_row] = evec[pivot_row], evec[k]
        end
        pivot = R[k, k]
        iszero(pivot) && throw(ArgumentError("APA Gram matrix is singular; increase regularization"))
        inv_pivot = one(T) / pivot
        for j in k:n
            R[k, j] *= inv_pivot
        end
        evec[k] *= inv_pivot
        for i in 1:n
            i == k && continue
            factor = R[i, k]
            iszero(factor) && continue
            for j in k:n
                R[i, j] -= factor * R[k, j]
            end
            evec[i] -= factor * evec[k]
        end
    end
    return evec
end

function _equalizer_slicer(y::T, constellation::AbstractVector{<:Number}) where {T<:Number}
    best_point = T(first(constellation))
    best_distance = abs2(y - best_point)
    for point in constellation
        candidate = T(point)
        distance = abs2(y - candidate)
        if distance < best_distance
            best_distance = distance
            best_point = candidate
        end
    end
    return best_point
end

function _equalizer_reference(step::Int, delay::Int, num_training::Int,
                               training_symbols::AbstractVector{<:Number}, y::T,
                               constellation::AbstractVector{<:Number}) where {T<:Number}
    if delay < step <= delay + num_training
        offset = step - delay - 1
        return T(training_symbols[firstindex(training_symbols) + offset])
    end
    return _equalizer_slicer(y, constellation)
end

_equalizer_plain_inner(a::Vector{T}, b::Vector{T}) where {T} = sum(a[k] * b[k] for k in eachindex(a, b))

function _equalize_kernel!(equalized::Vector{T}, err::Vector{T}, weights::Vector{T},
                            x::AbstractVector{<:Number}, training_symbols::AbstractVector{<:Number},
                            constellation::AbstractVector{<:Number}, reference_tap::Int, num_forward_taps::Int,
                            step_size::Real, ::Val{:lms}) where {T<:Number}
    num_taps = length(weights)
    num_feedback_taps = num_taps - num_forward_taps
    delay = reference_tap - 1
    num_training = length(training_symbols)
    mu = T(step_size)
    u = zeros(T, num_taps)
    forward_view = view(u, 1:num_forward_taps)
    feedback_view = view(u, (num_forward_taps + 1):num_taps)
    prev_ref = zero(T)

    for (step, sample) in enumerate(x)
        _equalizer_shift!(forward_view, T(sample))
        num_feedback_taps > 0 && _equalizer_shift!(feedback_view, prev_ref)
        y = _equalizer_output(weights, u)
        equalized[step] = y
        ref = _equalizer_reference(step, delay, num_training, training_symbols, y, constellation)
        e = ref - y
        err[step] = e
        for k in 1:num_taps
            weights[k] += mu * u[k] * conj(e)
        end
        prev_ref = ref
    end
    return equalized, err, weights
end

function _equalize_kernel!(equalized::Vector{T}, err::Vector{T}, weights::Vector{T},
                            x::AbstractVector{<:Number}, training_symbols::AbstractVector{<:Number},
                            constellation::AbstractVector{<:Number}, reference_tap::Int, num_forward_taps::Int,
                            forgetting_factor::Real, initial_inverse_correlation::Real,
                            ::Val{:rls}) where {T<:Number}
    num_taps = length(weights)
    num_feedback_taps = num_taps - num_forward_taps
    delay = reference_tap - 1
    num_training = length(training_symbols)
    lambda = T(forgetting_factor)
    lambda_inverse = one(T) / lambda
    u = zeros(T, num_taps)
    forward_view = view(u, 1:num_forward_taps)
    feedback_view = view(u, (num_forward_taps + 1):num_taps)
    prev_ref = zero(T)
    P = zeros(T, num_taps, num_taps)
    for k in 1:num_taps
        P[k, k] = T(initial_inverse_correlation)
    end
    UP = Vector{T}(undef, num_taps)
    Pu = Vector{T}(undef, num_taps)
    gain = Vector{T}(undef, num_taps)

    for (step, sample) in enumerate(x)
        _equalizer_shift!(forward_view, T(sample))
        num_feedback_taps > 0 && _equalizer_shift!(feedback_view, prev_ref)
        y = _equalizer_output(weights, u)
        equalized[step] = y
        ref = _equalizer_reference(step, delay, num_training, training_symbols, y, constellation)
        e = ref - y
        err[step] = e

        for j in 1:num_taps
            acc = zero(T)
            for k in 1:num_taps
                acc += conj(u[k]) * P[k, j]
            end
            UP[j] = acc
        end
        for i in 1:num_taps
            acc = zero(T)
            for k in 1:num_taps
                acc += P[i, k] * u[k]
            end
            Pu[i] = acc
        end

        denom = lambda + _equalizer_plain_inner(UP, u)
        for k in 1:num_taps
            gain[k] = Pu[k] / denom
            weights[k] += gain[k] * conj(e)
        end
        for i in 1:num_taps, j in 1:num_taps
            P[i, j] = (P[i, j] - gain[i] * UP[j]) * lambda_inverse
        end
        prev_ref = ref
    end
    return equalized, err, weights
end

function _equalize_kernel!(equalized::Vector{T}, err::Vector{T}, weights::Vector{T},
                            x::AbstractVector{<:Number}, training_symbols::AbstractVector{<:Number},
                            constellation::AbstractVector{<:Number}, reference_tap::Int, num_forward_taps::Int,
                            projection_order::Int, step_size::Real, regularization::Real,
                            ::Val{:apa}) where {T<:Number}
    num_taps = length(weights)
    num_feedback_taps = num_taps - num_forward_taps
    delay = reference_tap - 1
    num_training = length(training_symbols)
    mu = T(step_size)
    delta = T(regularization)
    N = projection_order
    u = zeros(T, num_taps)
    forward_view = view(u, 1:num_forward_taps)
    feedback_view = view(u, (num_forward_taps + 1):num_taps)
    prev_ref = zero(T)
    Uhist = zeros(T, num_taps, N)
    dhist = zeros(T, N)
    R = Matrix{T}(undef, N, N)
    evec = Vector{T}(undef, N)

    for (step, sample) in enumerate(x)
        _equalizer_shift!(forward_view, T(sample))
        num_feedback_taps > 0 && _equalizer_shift!(feedback_view, prev_ref)
        y = _equalizer_output(weights, u)
        equalized[step] = y
        ref = _equalizer_reference(step, delay, num_training, training_symbols, y, constellation)
        e = ref - y
        err[step] = e

        _equalizer_shift_matrix_columns!(Uhist, u)
        _equalizer_shift!(dhist, ref)

        for j in 1:N
            d_val = zero(T)
            for k in 1:num_taps
                d_val += conj(weights[k]) * Uhist[k, j]
            end
            evec[j] = conj(dhist[j] - d_val)
            for i in 1:N
                acc = zero(T)
                for k in 1:num_taps
                    acc += conj(Uhist[k, i]) * Uhist[k, j]
                end
                R[i, j] = acc
            end
            R[j, j] += delta
        end

        _equalizer_apa_solve!(R, evec)

        for k in 1:num_taps
            acc = zero(T)
            for j in 1:N
                acc += Uhist[k, j] * evec[j]
            end
            weights[k] += mu * acc
        end

        prev_ref = ref
    end
    return equalized, err, weights
end
