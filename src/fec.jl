const _FEC_MAX_CONSTRAINT_LENGTH = 24
const _FEC_MAX_OUTPUT_BITS = 16
const _FEC_MAX_SOFT_BITS = 52

"""
    poly2trellis(constraint_length, code_generator)

Build a rate-`1/n`, single-input-stream (`k=1`), non-recursive binary
convolutional-code trellis. `code_generator` is an `n`-length vector of
octal tap-connection integers, written as Julia octal literals (for example,
`[0o171, 0o133]` for the standard rate-1/2 `K=7` code).

Returns a plain `NamedTuple` (this library never uses structs for this):
`constraint_length`, `code_generator`, `k` (always `1`), `n`, `num_states`,
`num_input_symbols` (always `2`), `num_output_symbols`, `next_states`
(`num_states`x`2` matrix of 0-based next-state values), and `outputs`
(`num_states`x`2` matrix of decimal, MSB-first-packed output symbols).

`constraint_length` must be at least `2` and is capped at
`$_FEC_MAX_CONSTRAINT_LENGTH` to prevent accidental exponential allocation;
`num_states` doubles for each increment. A generator entry of `0` is legal.

Out of scope for this first version, deferred rather than assumed
unnecessary: multiple input streams (`k>1`), recursive/systematic-feedback
generators, and puncturing.
"""
function poly2trellis end

"""
    convenc(bits, trellis; initial_state=0)

Rate-`1/n` non-recursive convolutional encode. Returns `(code, final_state)`;
`code` has length `length(bits) * trellis.n`. This function never appends a
termination tail. For
[`vitdec`](@ref)'s `:term` mode (which assumes the encoder was driven back
to state `0`), the caller must append `trellis.constraint_length - 1` zero
bits to `bits` before calling. An empty `bits`
returns `(BitVector(), initial_state)` immediately.
"""
function convenc end

"""
    vitdec(code, trellis, traceback_depth; mode=:trunc, decision_type=:hard, num_soft_bits=1)

Viterbi decode. `mode` is `:trunc` (traceback starts from the state with the
lowest final path metric; ties favor the highest-index state) or `:term`
(traceback starts from state `0`, assuming a caller-appended zero tail).
This function never appends or strips that tail. `traceback_depth` is
validated for range but does not affect decoded bits in these full-block,
non-windowed modes. `:cont` mode
and puncturing/erasures are out of scope for this version. Returns
`length(code) / trellis.n * trellis.k` decoded bits always -- no
automatic tail-stripping for `:term`. An empty `code`
returns an empty `BitVector` immediately.

`decision_type` selects how `code` is interpreted. `:hard` expects entries
of exactly `0`/`1`. `:soft` expects `code` entries as integers in
`[0, 2^num_soft_bits - 1]`, where `0` is the most confident `0` and
`2^num_soft_bits - 1` is the most confident `1`. `:unquant` expects `code`
entries as arbitrary finite reals (e.g. raw matched-filter samples), where
`+1` represents logical `0` and `-1` logical `1`. `:llr` accepts arbitrary
finite log-likelihood ratios, with positive values favoring logical `0` and
negative values favoring logical `1`. Per-bit-position branch
metrics are computed once per received symbol into a
`trellis.num_output_symbols`-length table, then looked up during the
add-compare-select step -- for a candidate output symbol `s` and
bit-position `j` (`1`-indexed, matching this library's own MSB-first
`_bits_to_integer`/`_integer_to_bits!` convention: `code`'s `j`-th entry in
a symbol corresponds to bit `n-j` of `s`, counting from the LSB), a
candidate bit of `1` contributes `2^num_soft_bits - 1 - code[j]` (`:soft`)
or `(code[j]+1)^2` (`:unquant`); a candidate bit of `0` contributes
`code[j]` (`:soft`) or `(code[j]-1)^2` (`:unquant`). For `:llr`, a candidate
bit of `1` contributes `max(llr, 0)` and a candidate bit of `0` contributes
`max(-llr, 0)` after one block-global positive rescaling, which preserves all
path decisions while preventing finite extreme metrics from overflowing.
`:hard`'s metric
(`count_ones(xor(s, received_symbol))`, a whole-symbol Hamming distance)
is mathematically the same per-bit sum specialized to `0`/`1` inputs.
`num_soft_bits=1` with `decision_type=:soft` is therefore exactly
equivalent to `decision_type=:hard` on the same `0`/`1` data (`2^1-1=1`
collapses the soft formula to the hard one) -- a useful correctness
invariant.

`num_soft_bits` is validated unconditionally (`>= 1`, not `Bool`) and is
unused for `:hard` and `:unquant`. Unquantized inputs must be finite; the
scale-safe branch metric handles all finite floating-point magnitudes.
"""
function vitdec end

function poly2trellis(constraint_length::Integer, code_generator::AbstractVector{<:Integer})
    _validate_poly2trellis_arguments(constraint_length, code_generator)
    K = Int(constraint_length)
    generators = Int.(collect(code_generator))
    n = length(generators)
    num_states = 1 << (K - 1)
    next_states, outputs = _poly2trellis_build(K, generators)
    return (
        constraint_length=K,
        code_generator=generators,
        k=1,
        n=n,
        num_states=num_states,
        num_input_symbols=2,
        num_output_symbols=1 << n,
        next_states=next_states,
        outputs=outputs,
    )
end

function convenc(bits::AbstractVector{<:Integer}, trellis; initial_state::Integer=0)
    _validate_trellis_shape(trellis)
    _validate_convenc_arguments(trellis, initial_state)
    _validate_bit_vector(bits)

    isempty(bits) && return BitVector(), Int(initial_state)

    code = BitVector(undef, length(bits) * trellis.n)
    final_state = _convenc_kernel!(code, bits, trellis, Int(initial_state))
    return code, final_state
end

function vitdec(code::AbstractVector{<:Real}, trellis, traceback_depth::Integer; mode::Symbol=:trunc,
                 decision_type::Symbol=:hard, num_soft_bits::Integer=1)
    _validate_trellis_shape(trellis)
    length(code) % trellis.n == 0 || throw(ArgumentError("code length must be a multiple of trellis.n"))

    isempty(code) && return BitVector()

    symbol_count = length(code) ÷ trellis.n
    _validate_vitdec_arguments(traceback_depth, mode, decision_type, num_soft_bits, symbol_count)
    _validate_vitdec_code(code, decision_type, Int(num_soft_bits))

    decoded = BitVector(undef, symbol_count * trellis.k)
    _vitdec_kernel!(decoded, code, trellis, mode, decision_type, Int(num_soft_bits))
    return decoded
end

function _validate_poly2trellis_arguments(constraint_length::Integer, code_generator::AbstractVector{<:Integer})
    constraint_length isa Bool && throw(ArgumentError("constraint_length must be an integer, not Bool"))
    (constraint_length >= 2 && constraint_length <= _FEC_MAX_CONSTRAINT_LENGTH) ||
        throw(ArgumentError("constraint_length must be between 2 and $_FEC_MAX_CONSTRAINT_LENGTH"))

    isempty(code_generator) && throw(ArgumentError("code_generator must not be empty"))
    length(code_generator) <= _FEC_MAX_OUTPUT_BITS ||
        throw(ArgumentError("code_generator length must be <= $_FEC_MAX_OUTPUT_BITS"))

    limit = 1 << Int(constraint_length)
    for g in code_generator
        g isa Bool && throw(ArgumentError("code_generator entries must be integers, not Bool"))
        (g >= 0 && g < limit) || throw(ArgumentError("code_generator entries must be in [0, 2^constraint_length)"))
    end
    return nothing
end

function _poly2trellis_build(K::Int, generators::Vector{Int})
    mem = K - 1
    num_states = 1 << mem
    next_states = Matrix{Int}(undef, num_states, 2)
    outputs = Matrix{Int}(undef, num_states, 2)
    for state in 0:(num_states - 1)
        for u in 0:1
            next_states[state + 1, u + 1] = (state >> 1) | (u << (mem - 1))
            full_register = (u << mem) | state
            symbol = 0
            for g in generators
                symbol = (symbol << 1) | (count_ones(full_register & g) & 1)
            end
            outputs[state + 1, u + 1] = symbol
        end
    end
    return next_states, outputs
end

function _validate_trellis_shape(trellis)
    trellis.k == 1 || throw(ArgumentError("trellis.k must be 1"))
    trellis.n isa Integer && (1 <= trellis.n <= _FEC_MAX_OUTPUT_BITS) ||
        throw(ArgumentError("trellis.n must be in [1, $_FEC_MAX_OUTPUT_BITS]"))
    trellis.num_states isa Integer && (1 <= trellis.num_states <= typemax(Int32)) ||
        throw(ArgumentError("trellis.num_states must be in [1, typemax(Int32)]"))
    trellis.num_input_symbols == 2 || throw(ArgumentError("trellis.num_input_symbols must be 2"))
    trellis.num_output_symbols == (1 << Int(trellis.n)) ||
        throw(ArgumentError("trellis.num_output_symbols must equal 2^trellis.n"))
    size(trellis.next_states) == (trellis.num_states, 2) ||
        throw(ArgumentError("trellis.next_states must be num_states x 2"))
    size(trellis.outputs) == (trellis.num_states, 2) ||
        throw(ArgumentError("trellis.outputs must be num_states x 2"))
    all(state -> state isa Integer && 0 <= state < trellis.num_states, trellis.next_states) ||
        throw(ArgumentError("trellis.next_states entries must be in [0, num_states)"))
    all(output -> output isa Integer && 0 <= output < trellis.num_output_symbols, trellis.outputs) ||
        throw(ArgumentError("trellis.outputs entries must be in [0, num_output_symbols)"))
    return nothing
end

function _validate_convenc_arguments(trellis, initial_state::Integer)
    initial_state isa Bool && throw(ArgumentError("initial_state must be an integer, not Bool"))
    (initial_state >= 0 && initial_state < trellis.num_states) ||
        throw(ArgumentError("initial_state must be in [0, num_states)"))
    return nothing
end

function _convenc_kernel!(code::BitVector, bits::AbstractVector{<:Integer}, trellis, state::Int)
    output_start = firstindex(code)
    n = trellis.n
    step = 0
    for bit in bits
        input = Int(bit)
        output_symbol = trellis.outputs[state + 1, input + 1]
        _integer_to_bits!(code, output_start + step * n, n, output_symbol)
        state = trellis.next_states[state + 1, input + 1]
        step += 1
    end
    return state
end

function _validate_vitdec_arguments(traceback_depth::Integer, mode::Symbol, decision_type::Symbol,
                                     num_soft_bits::Integer, symbol_count::Int)
    mode in (:trunc, :term) || throw(ArgumentError("mode must be :trunc or :term (:cont is not implemented)"))
    traceback_depth isa Bool && throw(ArgumentError("traceback_depth must be an integer, not Bool"))
    (traceback_depth >= 1 && traceback_depth <= symbol_count) ||
        throw(ArgumentError("traceback_depth must be in [1, length(code) ÷ trellis.n]"))

    decision_type in (:hard, :soft, :unquant, :llr) ||
        throw(ArgumentError("decision_type must be :hard, :soft, :unquant, or :llr"))

    num_soft_bits isa Bool && throw(ArgumentError("num_soft_bits must be an integer, not Bool"))
    (1 <= num_soft_bits <= _FEC_MAX_SOFT_BITS) ||
        throw(ArgumentError("num_soft_bits must be in [1, $_FEC_MAX_SOFT_BITS]"))
    return nothing
end

function _validate_vitdec_code(code::AbstractVector{<:Real}, decision_type::Symbol, num_soft_bits::Int)
    if decision_type === :hard
        for value in code
            (value == 0 || value == 1) ||
                throw(ArgumentError("code must contain only 0 or 1 for decision_type=:hard"))
        end
    elseif decision_type === :soft
        limit = (1 << num_soft_bits) - 1
        for value in code
            (isfinite(value) && value == floor(value) && value >= 0 && value <= limit) ||
                throw(ArgumentError("code must contain integers in [0, 2^num_soft_bits - 1] for decision_type=:soft"))
        end
    elseif decision_type === :unquant
        for value in code
            (isfinite(value) && isfinite(Float64(value))) ||
                throw(ArgumentError("code must be finite and representable as Float64 for decision_type=:unquant"))
        end
    else
        for value in code
            (isfinite(value) && isfinite(Float64(value))) ||
                throw(ArgumentError("code must be finite and representable as Float64 for decision_type=:llr"))
        end
    end
    return nothing
end

function _vitdec_branch_metric_table!(table::Vector{Float64}, code::AbstractVector{<:Real}, start::Int, n::Int,
                                       num_output_symbols::Int, decision_type::Symbol, num_soft_bits::Int,
                                       llr_scale::Float64)
    if decision_type === :hard
        received = 0
        for j in 1:n
            received = (received << 1) | Int(code[start + j - 1])
        end
        @inbounds for s in 0:(num_output_symbols - 1)
            table[s + 1] = count_ones(xor(s, received))
        end
    elseif decision_type === :soft
        max_soft = Float64((1 << num_soft_bits) - 1)
        @inbounds for s in 0:(num_output_symbols - 1)
            acc = 0.0
            for j in 1:n
                bit = (s >> (n - j)) & 1
                cw = Float64(code[start + j - 1])
                acc += bit == 1 ? (max_soft - cw) : cw
            end
            table[s + 1] = acc
        end
    elseif decision_type === :unquant
        max_abs = 0.0
        for j in 1:n
            max_abs = max(max_abs, abs(Float64(code[start + j - 1])))
        end
        if max_abs <= sqrt(floatmax(Float64) / n)
            @inbounds for s in 0:(num_output_symbols - 1)
                acc = 0.0
                for j in 1:n
                    bit = (s >> (n - j)) & 1
                    cw = Float64(code[start + j - 1])
                    acc += bit == 1 ? (cw + 1)^2 : (cw - 1)^2
                end
                table[s + 1] = acc
            end
        else
            inverse_scale = inv(max_abs)
            @inbounds for s in 0:(num_output_symbols - 1)
                acc = 0.0
                for j in 1:n
                    bit = (s >> (n - j)) & 1
                    cw = Float64(code[start + j - 1]) * inverse_scale
                    acc -= bit == 1 ? -cw : cw
                end
                table[s + 1] = acc
            end
        end
    else
        @inbounds for s in 0:(num_output_symbols - 1)
            acc = 0.0
            for j in 1:n
                bit = (s >> (n - j)) & 1
                llr = Float64(code[start + j - 1]) / llr_scale
                acc += bit == 1 ? max(llr, 0.0) : max(-llr, 0.0)
            end
            table[s + 1] = acc
        end
    end
    return table
end

function _vitdec_kernel!(decoded::BitVector, code::AbstractVector{<:Real}, trellis, mode::Symbol,
                          decision_type::Symbol, num_soft_bits::Int)
    num_states = trellis.num_states
    n = trellis.n
    num_output_symbols = trellis.num_output_symbols
    next_states = trellis.next_states
    outputs = trellis.outputs
    input_start = firstindex(code)
    symbol_count = length(code) ÷ n
    sentinel = Inf

    metric = fill(sentinel, num_states)
    metric[1] = 0.0
    new_metric = Vector{Float64}(undef, num_states)
    predecessor_state = Matrix{Int32}(undef, num_states, symbol_count)
    predecessor_input = falses(num_states, symbol_count)
    branch_metric_table = Vector{Float64}(undef, num_output_symbols)
    llr_scale = decision_type === :llr ? max(1.0, maximum(abs ∘ Float64, code)) : 1.0

    for step in 1:symbol_count
        start = input_start + (step - 1) * n
        _vitdec_branch_metric_table!(branch_metric_table, code, start, n, num_output_symbols, decision_type,
                                      num_soft_bits, llr_scale)
        fill!(new_metric, sentinel)
        @inbounds for state in 0:(num_states - 1)
            current = metric[state + 1]
            !isfinite(current) && continue
            for u in 0:1
                ns = next_states[state + 1, u + 1]
                candidate_output = outputs[state + 1, u + 1]
                branch_metric = branch_metric_table[candidate_output + 1]
                candidate = current + branch_metric
                if candidate < new_metric[ns + 1]
                    new_metric[ns + 1] = candidate
                    predecessor_state[ns + 1, step] = state
                    predecessor_input[ns + 1, step] = isodd(u)
                end
            end
        end
        metric, new_metric = new_metric, metric
    end

    final_state = 0
    if mode === :trunc
        best_metric = metric[1]
        for state in 1:(num_states - 1)
            if metric[state + 1] <= best_metric
                best_metric = metric[state + 1]
                final_state = state
            end
        end
    end

    isfinite(metric[final_state + 1]) || throw(ArgumentError("Viterbi decoder found no finite survivor path"))

    output_start = firstindex(decoded)
    state = final_state
    @inbounds for step in symbol_count:-1:1
        decoded[output_start + step - 1] = predecessor_input[state + 1, step]
        state = Int(predecessor_state[state + 1, step])
    end
    return decoded
end
