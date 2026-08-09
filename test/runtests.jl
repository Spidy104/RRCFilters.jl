using RRCFilters
using Random
using Statistics
using Test

struct ZeroBasedVector{T} <: AbstractVector{T}
    data::Vector{T}
end

Base.size(vector::ZeroBasedVector) = (length(vector.data),)
Base.axes(vector::ZeroBasedVector) = (0:(length(vector.data) - 1),)
Base.IndexStyle(::Type{<:ZeroBasedVector}) = IndexCartesian()
Base.getindex(vector::ZeroBasedVector, index::Int) = (checkbounds(vector, index); vector.data[index + 1])

function reference_filter(beta, span, sps, shape)
    setprecision(BigFloat, 256) do
        b = BigFloat(beta)
        count = span * sps
        half_count = count ÷ 2
        h = Vector{BigFloat}(undef, count + 1)

        for (index, n) in enumerate(-half_count:half_count)
            t = BigFloat(n) / BigFloat(sps)
            if iszero(b)
                h[index] = sinc(t)
            elseif shape === :normal
                scaled_time = 2b * t
                h[index] = iszero(t) ? one(BigFloat) :
                    isapprox(abs(scaled_time), one(BigFloat); rtol=8eps(BigFloat), atol=8eps(BigFloat)) ?
                    b / 2 * sin(big(pi) / (2b)) :
                    sinc(t) * cos(big(pi) * b * t) / (1 - scaled_time^2)
            else
                scaled_time = 4b * t
                if iszero(t)
                    h[index] = 1 - b + 4b / big(pi)
                elseif isapprox(abs(scaled_time), one(BigFloat); rtol=8eps(BigFloat), atol=8eps(BigFloat))
                    angle = big(pi) / (4b)
                    h[index] = b / sqrt(big(2)) * ((1 + 2 / big(pi)) * sin(angle) + (1 - 2 / big(pi)) * cos(angle))
                else
                    h[index] = (sin(big(pi) * t * (1 - b)) + 4b * t * cos(big(pi) * t * (1 + b))) /
                               (big(pi) * t * (1 - scaled_time^2))
                end
            end
        end

        h ./= sqrt(sum(abs2, h))
        return h
    end
end

function reference_upfirdn(x::AbstractVector{<:Number}, h::AbstractVector{<:Number}, p::Integer, q::Integer)
    nx = length(x)
    nh = length(h)
    isempty(x) && return Complex{BigFloat}[]

    upsampled = zeros(Complex{BigFloat}, (nx - 1) * p + 1)
    for (index, sample) in enumerate(x)
        upsampled[(index - 1) * p + 1] = Complex{BigFloat}(sample)
    end

    full = zeros(Complex{BigFloat}, length(upsampled) + nh - 1)
    for (hi, tap) in enumerate(h), (ui, sample) in enumerate(upsampled)
        full[ui + hi - 1] += Complex{BigFloat}(tap) * sample
    end

    output_length = cld((nx - 1) * p + nh, q)
    return [full[(k - 1) * q + 1] for k in 1:output_length]
end

const RRC_REFERENCE = [
    -0.01877334400452212, 0.003013558666422048, 0.03267723454625457,
    0.047093583340252425, 0.02654951770229082, -0.027522224030017674,
    -0.08522487501708324, -0.09944743599192618, -0.032147267398314666,
    0.11903714817925416, 0.31117641155547854, 0.47200389556543276,
    0.5346320701781558, 0.47200389556543276, 0.31117641155547854,
    0.11903714817925416, -0.032147267398314666, -0.09944743599192618,
    -0.08522487501708324, -0.027522224030017674, 0.02654951770229082,
    0.047093583340252425, 0.03267723454625457, 0.003013558666422048,
    -0.01877334400452212,
]

const RC_REFERENCE = [
    0.0, 0.026393742324270428, 0.04477975389257778, 0.03798165105400917,
    -1.5827244461189454e-17, -0.055344691535841976, -0.09595661548409526,
    -0.08486603301188494, 0.0, 0.15010390207108237, 0.32432366750157204,
    0.46372758906276906, 0.5169570351944414, 0.46372758906276906,
    0.32432366750157204, 0.15010390207108237, 0.0, -0.08486603301188494,
    -0.09595661548409526, -0.055344691535841976, -1.5827244461189454e-17,
    0.03798165105400917, 0.04477975389257778, 0.026393742324270428, 0.0,
]

@testset "RRCFilters" begin
    @testset "input validation" begin
        @test_throws ArgumentError rcosdesign(NaN, 6, 4)
        @test_throws ArgumentError rcosdesign(Inf, 6, 4)
        @test_throws ArgumentError rcosdesign(-0.1, 6, 4)
        @test_throws ArgumentError rcosdesign(1.1, 6, 4)
        @test_throws ArgumentError rcosdesign(true, 6, 4)
        @test_throws ArgumentError rcosdesign(0.25, true, 4)
        @test_throws ArgumentError rcosdesign(0.25, 6, false)
        @test_throws ArgumentError rcosdesign(0.25, 0, 4)
        @test_throws ArgumentError rcosdesign(0.25, 6, 0)
        @test_throws ArgumentError rcosdesign(0.25, -6, 4)
        @test_throws ArgumentError rcosdesign(0.25, 6, -4)
        @test_throws ArgumentError rcosdesign(0.25, 3, 3)
        @test_throws ArgumentError rcosdesign(0.25, typemax(Int), 2)
        @test_throws ArgumentError rcosdesign(0.25, 6, 4; shape=:invalid)
        @test_throws MethodError rcosdesign(0.25, 6.0, 4)
        @test_throws MethodError rcosdesign(0.25, 6, 4.0)
        @test_throws TypeError rcosdesign(0.25, 6, 4; shape="sqrt")
    end

    @testset "fixed reference designs" begin
        @test rcosdesign(0.25, 6, 4; shape=:sqrt) ≈ RRC_REFERENCE rtol=2e-15
        @test rcosdesign(0.25, 6, 4; shape=:normal) ≈ RC_REFERENCE rtol=2e-15 atol=2e-15
    end

    @testset "filter invariants" begin
        for beta in (0.0, 0.1, 0.25, 0.5, 0.75, 1.0), span in (2, 4, 6, 8), sps in (2, 4, 8), shape in (:sqrt, :normal)
            h = rcosdesign(beta, span, sps; shape)
            @test length(h) == span * sps + 1
            @test h == reverse(h)
            @test all(isfinite, h)
            @test sum(abs2, h) ≈ 1.0 rtol=16eps(Float64)
        end
    end

    @testset "high precision reference" begin
        for beta in (0.1, 0.25, 0.35, 0.5, 1.0), shape in (:sqrt, :normal)
            expected = reference_filter(beta, 6, 4, shape)
            @test rcosdesign(beta, 6, 4; shape) ≈ Float64.(expected) rtol=2e-14 atol=2e-14
        end
    end

    @testset "special cases" begin
        rc_zero = rcosdesign(0.0, 6, 4; shape=:normal)
        rrc_zero = rcosdesign(0.0, 6, 4; shape=:sqrt)
        @test rc_zero ≈ rrc_zero rtol=2e-15
        @test maximum(abs, rc_zero[[1, 5, 9, 17, 21, 25]]) < 1e-14
        @test rcosdesign(0.5, 4, 4; shape=:normal) ≈ Float64.(reference_filter(0.5, 4, 4, :normal)) rtol=2e-14
        @test rcosdesign(0.25, 4, 4; shape=:sqrt) ≈ Float64.(reference_filter(0.25, 4, 4, :sqrt)) rtol=2e-14
    end

    @testset "near singularities" begin
        # These roll-offs place a sampled point 1e-12 from each removable
        # singularity, where the direct formulas suffer cancellation.
        rc_beta = 1 / 3 + 1e-12
        rrc_beta = 0.25 + 1e-12
        @test rcosdesign(rc_beta, 6, 2; shape=:normal) ≈ Float64.(reference_filter(rc_beta, 6, 2, :normal)) rtol=1e-10 atol=1e-12
        @test rcosdesign(rrc_beta, 6, 4; shape=:sqrt) ≈ Float64.(reference_filter(rrc_beta, 6, 4, :sqrt)) rtol=1e-10 atol=1e-12
    end

    @testset "numeric types and compiler behavior" begin
        h32 = rcosdesign(Float32(0.25), 6, 4)
        @test eltype(h32) === Float32
        @test sum(abs2, h32) ≈ Float32(1) rtol=32eps(Float32)
        hbig = rcosdesign(big"0.25", 6, 4)
        @test eltype(hbig) === BigFloat
        @test sum(abs2, hbig) ≈ big"1.0" rtol=32eps(BigFloat)
        @test @inferred(rcosdesign(0.25, 6, 4)) isa Vector{Float64}
        @test isempty(Test.detect_ambiguities(RRCFilters))
    end

    @testset "awgn" begin
        @testset "real and complex SNR" begin
            sample_count = 100_000
            real_signal = ones(sample_count)
            real_noisy = awgn(Xoshiro(1), real_signal, 10.0)
            real_noise = real_noisy - real_signal
            @test 10 * log10(mean(abs2, real_signal) / mean(abs2, real_noise)) ≈ 10.0 atol=0.1

            negative_snr = awgn(Xoshiro(2), real_signal, -5.0)
            negative_noise = negative_snr - real_signal
            @test 10 * log10(mean(abs2, real_signal) / mean(abs2, negative_noise)) ≈ -5.0 atol=0.1

            complex_signal = fill((1 + im) / sqrt(2), sample_count)
            complex_noisy = awgn(Xoshiro(3), complex_signal, 8.0)
            complex_noise = complex_noisy - complex_signal
            @test 10 * log10(mean(abs2, complex_signal) / mean(abs2, complex_noise)) ≈ 8.0 atol=0.1
            @test mean(abs2, real.(complex_noise)) ≈ mean(abs2, imag.(complex_noise)) rtol=0.03
        end

        @testset "RNG and signal power" begin
            signal = fill(2.0, 10_000)
            override = awgn(Xoshiro(4), signal, 10.0; signal_power=1.0) - signal
            measured = awgn(Xoshiro(4), signal, 10.0) - signal
            @test measured ≈ 2 .* override
            @test awgn(Xoshiro(5), signal, 3.0) == awgn(Xoshiro(5), signal, 3.0)
            advancing_rng = Xoshiro(6)
            @test awgn(advancing_rng, signal, 3.0) != awgn(advancing_rng, signal, 3.0)
            @test awgn(Xoshiro(7), zeros(Float32, 4), 10.0) == zeros(Float32, 4)
            @test awgn(Xoshiro(8), signal, 10.0; signal_power=0.0) == signal
            @test length(awgn([1.0], 10.0)) == 1
        end

        @testset "types and immutability" begin
            real32 = Float32[1, -1]
            real16 = Float16[1, -1]
            real64 = Float64[1, -1]
            complex32 = ComplexF32[1 + im, -1 - im]
            complex64 = ComplexF64[1 + im, -1 - im]
            integers = Int[1, -1]
            complex_integers = Complex{Int}[1 + im, -1 - im]
            rationals = Rational{Int}[1 // 2, -1 // 2]
            bigfloats = BigFloat[big"1.0", big"-1.0"]
            abstract_reals = Number[1, -1]
            abstract_complex = Complex[1 + im, -1 - im]
            mixed_numbers = Number[1, 1 + im]
            abstract_bigfloats = AbstractFloat[big"1.0", big"-1.0"]
            mixed_big_numbers = Number[big"1.0", 1 + im]
            original = copy(real64)

            @test eltype(awgn(Xoshiro(9), real32, 10.0)) === Float32
            @test eltype(awgn(Xoshiro(9), real16, 10.0)) === Float16
            @test eltype(awgn(Xoshiro(10), real64, 10.0)) === Float64
            @test eltype(awgn(Xoshiro(11), complex32, 10.0)) === ComplexF32
            @test eltype(awgn(Xoshiro(12), complex64, 10.0)) === ComplexF64
            @test eltype(awgn(Xoshiro(13), integers, 10.0)) === Float64
            @test eltype(awgn(Xoshiro(14), complex_integers, 10.0)) === ComplexF64
            @test eltype(awgn(Xoshiro(15), rationals, 10.0)) === Float64
            @test eltype(awgn(Xoshiro(16), bigfloats, 10.0)) === BigFloat
            @test eltype(awgn(Xoshiro(17), abstract_reals, 10.0)) === Float64
            @test eltype(awgn(Xoshiro(18), abstract_complex, 10.0)) === ComplexF64
            @test eltype(awgn(Xoshiro(19), mixed_numbers, 10.0)) === ComplexF64
            @test eltype(awgn(Xoshiro(19), abstract_bigfloats, 10.0)) === BigFloat
            @test eltype(awgn(Xoshiro(19), mixed_big_numbers, 10.0)) === Complex{BigFloat}
            @test all(isfinite, awgn(Xoshiro(20), [typemax(Int)], 0.0))
            @test length(awgn(Xoshiro(15), real64, 10.0)) == length(real64)
            awgn(Xoshiro(21), real64, 10.0)
            @test real64 == original
            @test @inferred(awgn(Xoshiro(22), real64, 10.0)) isa Vector{Float64}
        end

        @testset "extreme numeric ranges" begin
            tiny_signal = [1e-200]
            @test awgn(Xoshiro(23), tiny_signal, 0.0) != tiny_signal

            finite_mean_signal = fill(sqrt(0.75 * floatmax(Float64)), 2)
            @test all(isfinite, awgn(Xoshiro(24), finite_mean_signal, 0.0))
            @test all(isfinite, awgn(Xoshiro(25), ones(Float16, 100_000), Float16(10)))

            @test !iszero(only(awgn(Xoshiro(26), [0.0], 4000.0; signal_power=1.0)))
            @test isfinite(only(awgn(Xoshiro(27), [0.0], -4000.0; signal_power=1.0)))
            @test isfinite(only(awgn(Xoshiro(28), [0.0], 1000.0; signal_power=big"1e400")))

            @test_throws ArgumentError awgn(Xoshiro(29), [big(10)^400], 10.0)
            @test_throws ArgumentError awgn(Xoshiro(30), [1 // big(10)^1000], 10.0)
        end

        @testset "validation" begin
            @test_throws ArgumentError awgn(Xoshiro(18), Float64[], 10.0)
            @test_throws ArgumentError awgn(Xoshiro(18), [NaN], 10.0)
            @test_throws ArgumentError awgn(Xoshiro(18), [Inf], 10.0)
            @test_throws ArgumentError awgn(Xoshiro(18), [-Inf], 10.0)
            @test_throws ArgumentError awgn(Xoshiro(18), ComplexF64[1 + Inf * im], 10.0)
            @test_throws ArgumentError awgn(Xoshiro(18), [1.0], NaN)
            @test_throws ArgumentError awgn(Xoshiro(18), [1.0], Inf)
            @test_throws ArgumentError awgn(Xoshiro(18), [1.0], -Inf)
            @test_throws ArgumentError awgn(Xoshiro(18), [1.0], true)
            @test_throws ArgumentError awgn(Xoshiro(18), [1.0], 10.0; signal_power=NaN)
            @test_throws ArgumentError awgn(Xoshiro(18), [1.0], 10.0; signal_power=Inf)
            @test_throws ArgumentError awgn(Xoshiro(18), [1.0], 10.0; signal_power=-Inf)
            @test_throws ArgumentError awgn(Xoshiro(18), [1.0], 10.0; signal_power=-1.0)
            @test_throws ArgumentError awgn(Xoshiro(18), [1.0], 10.0; signal_power=true)
            @test_throws ArgumentError awgn(Xoshiro(18), [1.0], -10_000.0)
            @test all(isfinite, awgn(Xoshiro(18), [floatmax(Float64)], 0.0))
            @test_throws MethodError awgn(Xoshiro(18), ["not a signal"], 10.0)
        end
    end

    @testset "metric foundations" begin
        @test RRCFilters._metric_bits_per_symbol(2) == 1
        @test RRCFilters._metric_bits_per_symbol(4) == 2
        @test RRCFilters._metric_bits_per_symbol(16) == 4
        @test RRCFilters._metric_bits_per_symbol(256) == 8
        @test_throws ArgumentError RRCFilters._metric_bits_per_symbol(true)
        @test_throws ArgumentError RRCFilters._metric_bits_per_symbol(1)
        @test_throws ArgumentError RRCFilters._metric_bits_per_symbol(12)
        @test_throws ArgumentError RRCFilters._metric_bits_per_symbol(4.0)

        metric_bits = BitVector([false, true, true, false])
        zero_based_metric_bits = ZeroBasedVector(collect(metric_bits))
        @test RRCFilters._validate_metric_vector_lengths(metric_bits, zero_based_metric_bits, "tx_bits", "rx_bits") == 4
        @test_throws ArgumentError RRCFilters._validate_metric_vector_lengths(metric_bits, BitVector([false]), "tx_bits", "rx_bits")
        @test_throws ArgumentError RRCFilters._validate_metric_vector_lengths(BitVector(), BitVector(), "tx_bits", "rx_bits")

        @test RRCFilters._validate_biterr_configuration(:bit, nothing, 4) == 1
        @test RRCFilters._validate_biterr_configuration(:symbol, 16, 8) == 4
        @test_throws ArgumentError RRCFilters._validate_biterr_configuration(:bit, 4, 4)
        @test_throws ArgumentError RRCFilters._validate_biterr_configuration(:symbol, nothing, 4)
        @test_throws ArgumentError RRCFilters._validate_biterr_configuration(:symbol, 16, 6)
        @test_throws ArgumentError RRCFilters._validate_biterr_configuration(:invalid, nothing, 4)
        @test_throws ArgumentError RRCFilters._validate_biterr_configuration("bit", nothing, 4)

        @test RRCFilters._validate_evm_mode(:percent) === :percent
        @test RRCFilters._validate_evm_mode(:dB) === :dB
        @test_throws ArgumentError RRCFilters._validate_evm_mode(:rms)
        @test_throws ArgumentError RRCFilters._validate_evm_mode("percent")
    end

    @testset "biterr and evm" begin
        tx = BitVector([false, false, true, true])
        rx = BitVector([false, true, false, true])
        @test biterr(tx, tx) == (0, 0.0)
        @test biterr(tx, rx) == (2, 0.5)
        @test biterr(tx, .!tx) == (4, 1.0)
        @test biterr(tx, rx; mode=:symbol, M=4) == (2, 1.0)
        for M in (2, 4, 16, 64, 256)
            bits_per_symbol = trailing_zeros(M)
            grouped_tx = falses(2 * bits_per_symbol)
            grouped_rx = copy(grouped_tx)
            grouped_rx[1] = true
            grouped_rx[bits_per_symbol + 1] = true
            @test biterr(grouped_tx, grouped_rx) == (2, 1 / bits_per_symbol)
            @test biterr(grouped_tx, grouped_rx; mode=:symbol, M) == (2, 1.0)
        end
        @test_throws ArgumentError biterr(tx, BitVector([false]))
        @test_throws ArgumentError biterr(BitVector(), BitVector())
        @test_throws ArgumentError biterr(tx, rx; mode=:symbol)
        @test_throws ArgumentError biterr(tx, rx; mode=:bit, M=4)
        @test_throws ArgumentError biterr(tx, rx; mode=:symbol, M=true)
        @test_throws ArgumentError biterr(tx, rx; mode=:symbol, M=1)
        @test_throws ArgumentError biterr(tx, rx; mode=:symbol, M=8)
        @test_throws ArgumentError biterr(tx, rx; mode=:invalid)
        @test_throws ArgumentError biterr(tx, rx; mode="bit")

        zero_based_tx = ZeroBasedVector(Bool[false, false, true, true])
        zero_based_rx = ZeroBasedVector(Bool[false, true, false, true])
        @test biterr(zero_based_tx, zero_based_rx) == (2, 0.5)
        @test biterr(@view(tx[:]), @view(rx[:]); mode=:symbol, M=4) == (2, 1.0)
        @test @inferred(biterr(tx, rx)) == (2, 0.5)

        ideal = ComplexF64[1 + im, -1 - im]
        received = ComplexF64[1.1 + 0.9im, -0.8 - 1.2im]
        @test evm(ideal, received) > 0.0
        @test evm(ComplexF64[1 + 0im], ComplexF64[-1 + 0im]) == 200.0
        @test_throws ArgumentError evm(ideal, ComplexF64[1 + 0im])
        @test_throws ArgumentError evm(ComplexF64[], ComplexF64[])
        @test_throws ArgumentError evm(ComplexF64[0 + 0im], ComplexF64[0 + 0im])
        @test_throws ArgumentError evm(ComplexF64[NaN + 0im], ComplexF64[0 + 0im])
        @test_throws ArgumentError evm(ComplexF64[1 + 0im], ComplexF64[Inf + 0im])
        @test_throws ArgumentError evm(ideal, received; mode=:rms)
        @test_throws ArgumentError evm(ideal, received; mode="percent")

        zero_based_ideal = ZeroBasedVector(ComplexF64[1 + 0im, -1 + 0im])
        zero_based_received = ZeroBasedVector(ComplexF64[1.1 + 0im, -0.9 + 0im])
        @test evm(zero_based_ideal, zero_based_received) ≈ 10.0 rtol=8eps(Float64)
        @test evm(@view(ideal[:]), @view(received[:])) == evm(ideal, received)
        ideal32 = ComplexF32[1 + 0im, -1 + 0im]
        received32 = ComplexF32[1.1f0 + 0im, -0.9f0 + 0im]
        @test evm(ideal32, received32) ≈ 10.0 rtol=16eps(Float32)
        extreme_ideal = ComplexF64[floatmax(Float64) / 4 + 0im, -floatmax(Float64) / 4 + 0im]
        extreme_received = ComplexF64[0.9 * extreme_ideal[1], 0.9 * extreme_ideal[2]]
        @test evm(extreme_ideal, extreme_received) ≈ 10.0 rtol=32eps(Float64)

        largest = floatmax(Float64)
        smallest = nextfloat(0.0)
        opposite_ideal = ComplexF64[complex(largest, 0.0)]
        opposite_received = ComplexF64[complex(-largest, 0.0)]
        @test evm(opposite_ideal, opposite_received) == 200.0
        @test evm(opposite_ideal, opposite_received; mode=:dB) ≈ 20 * log10(2.0) rtol=16eps(Float64)

        wide_ideal = ComplexF64[complex(smallest, 0.0)]
        wide_received = ComplexF64[complex(largest, 0.0)]
        @test_throws ArgumentError evm(wide_ideal, wide_received)
        @test evm(wide_ideal, wide_received; mode=:dB) ≈ 12_631.218618060651 rtol=16eps(Float64)

        tiny_error_ideal = ComplexF64[complex(largest, 0.0)]
        tiny_error_received = ComplexF64[complex(largest, smallest)]
        @test_throws ArgumentError evm(tiny_error_ideal, tiny_error_received)
        @test evm(tiny_error_ideal, tiny_error_received; mode=:dB) ≈ -12_631.218618060651 rtol=16eps(Float64)

        huge_integer = big(10)^400
        bigint_ideal = Complex{BigInt}[complex(huge_integer, big(0))]
        bigint_received = Complex{BigInt}[complex(2 * huge_integer, big(0))]
        @test evm(bigint_ideal, bigint_received) ≈ big"100.0" rtol=32eps(BigFloat)
        @test evm(bigint_ideal, bigint_received; mode=:dB) == big"0.0"

        bigfloat_ideal = Complex{BigFloat}[complex(big"1e1000", big"0.0"), complex(-big"1e1000", big"0.0")]
        bigfloat_received = bigfloat_ideal .* big"1.1"
        @test evm(bigfloat_ideal, bigfloat_received) ≈ big"10.0" rtol=32eps(BigFloat)
        @test evm(bigfloat_ideal, bigfloat_received) isa BigFloat
        @test @inferred(evm(ideal, received)) == evm(ideal, received)
    end

    @testset "QAM link metrics" begin
        function balanced_qam_bits(M, symbol_count)
            parameters = RRCFilters._qam_parameters(M)
            bits = falses(symbol_count * parameters.bits_per_symbol)
            for symbol_index in 1:symbol_count
                label = (symbol_index - 1) % M
                RRCFilters._integer_to_bits!(bits, (symbol_index - 1) * parameters.bits_per_symbol + 1, parameters.bits_per_symbol, label)
            end
            return bits
        end

        for (M, low_snr, high_snr, seed) in ((16, 12.0, 24.0, 101), (64, 16.0, 28.0, 102), (256, 20.0, 32.0, 103))
            bits = balanced_qam_bits(M, 16_384)
            symbols = qammod(bits, M)
            low_noisy = awgn(Xoshiro(seed), symbols, low_snr)
            high_noisy = awgn(Xoshiro(seed), symbols, high_snr)
            low_bits = qamdemod(low_noisy, M)
            high_bits = qamdemod(high_noisy, M)
            low_errors, low_ber = biterr(bits, low_bits)
            high_errors, high_ber = biterr(bits, high_bits)
            low_symbol_errors, low_ser = biterr(bits, low_bits; mode=:symbol, M)
            high_symbol_errors, high_ser = biterr(bits, high_bits; mode=:symbol, M)
            low_evm = evm(symbols, low_noisy)
            high_evm = evm(symbols, high_noisy)

            @test mean(abs2, symbols) ≈ 1.0 rtol=8eps(Float64)
            @test low_evm ≈ 100 * 10^(-low_snr / 20) rtol=0.06
            @test high_evm ≈ 100 * 10^(-high_snr / 20) rtol=0.06
            @test evm(symbols, low_noisy; mode=:dB) ≈ 20 * log10(low_evm / 100) rtol=16eps(Float64)
            @test high_evm < low_evm
            @test high_ber <= low_ber
            @test high_ser <= low_ser
            @test low_errors >= high_errors
            @test low_symbol_errors >= high_symbol_errors
            @test 0.0 <= low_ber <= 1.0
            @test 0.0 <= low_ser <= 1.0
        end
    end

    @testset "QAM helper foundations" begin
        @test RRCFilters._qam_parameters(4) == (M=4, bits_per_symbol=2, kind=:square, axis_bits=1, side_length=2, scale=sqrt(0.5))
        @test RRCFilters._qam_parameters(16).bits_per_symbol == 4
        @test RRCFilters._qam_parameters(64).axis_bits == 3
        @test RRCFilters._qam_parameters(256).side_length == 16
        @test RRCFilters._qam_parameters(4).kind === :square
        @test RRCFilters._qam_parameters(8).kind === :nonsquare
        @test RRCFilters._qam_parameters(32).kind === :nonsquare
        @test_throws ArgumentError RRCFilters._qam_parameters(true)
        @test_throws ArgumentError RRCFilters._qam_parameters(2)
        @test_throws ArgumentError RRCFilters._qam_parameters(12)
        @test_throws ArgumentError RRCFilters._qam_parameters(2048)

        for value in 0:255
            @test RRCFilters._gray_to_binary(RRCFilters._binary_to_gray(value)) == value
        end

        bits = BitVector([false, true, true, false])
        @test RRCFilters._bits_to_integer(bits, 1, 4) == 6
        recovered = falses(4)
        RRCFilters._integer_to_bits!(recovered, 1, 4, 6)
        @test recovered == bits
        @test RRCFilters._qam_axis_level(0, 4) == -3
        @test RRCFilters._qam_axis_level(3, 4) == 3
        @test RRCFilters._qam_nearest_axis_position(-2.0, 4) == 1
        @test RRCFilters._qam_nearest_axis_position(2.0, 4) == 3
        @test_throws ArgumentError RRCFilters._validate_bit_vector([0, 2])

    end

    @testset "QAM modulation and demodulation" begin
        function qam_label_bits(M)
            parameters = RRCFilters._qam_parameters(M)
            bits = falses(M * parameters.bits_per_symbol)
            for label in 0:(M - 1)
                RRCFilters._integer_to_bits!(bits, label * parameters.bits_per_symbol + 1, parameters.bits_per_symbol, label)
            end
            return bits
        end

        for M in (4, 8, 16, 32, 64, 128, 256, 512), gray in (true, false)
            parameters = RRCFilters._qam_parameters(M)
            bits = qam_label_bits(M)
            original_bits = copy(bits)
            symbols = qammod(bits, M; gray)

            @test length(symbols) == M
            @test length(unique(symbols)) == M
            @test mean(symbols) ≈ 0.0 + 0.0im atol=4eps(Float64)
            @test mean(abs2, symbols) ≈ 1.0 rtol=8eps(Float64)
            @test qamdemod(symbols, M; gray) == bits
            @test qammod(@view(bits[:]), M; gray) == symbols
            @test qamdemod(@view(symbols[:]), M; gray) == bits
            @test bits == original_bits

            reference = qammod(bits[1:parameters.bits_per_symbol], M; gray)
            perturbed = reference .+ 0.1 * parameters.scale * (0.25 + 0.25im)
            @test qamdemod(perturbed, M; gray) == bits[1:parameters.bits_per_symbol]
        end

        for M in (4, 16, 64, 256)
            parameters = RRCFilters._qam_parameters(M)
            for i_position in 0:(parameters.side_length - 2), q_position in 0:(parameters.side_length - 1)
                left = parameters.scale * complex(RRCFilters._qam_axis_level(i_position, parameters.side_length), RRCFilters._qam_axis_level(q_position, parameters.side_length))
                right = parameters.scale * complex(RRCFilters._qam_axis_level(i_position + 1, parameters.side_length), RRCFilters._qam_axis_level(q_position, parameters.side_length))
                @test sum(qamdemod([left], M; gray=true) .!= qamdemod([right], M; gray=true)) == 1
            end
            for i_position in 0:(parameters.side_length - 1), q_position in 0:(parameters.side_length - 2)
                lower = parameters.scale * complex(RRCFilters._qam_axis_level(i_position, parameters.side_length), RRCFilters._qam_axis_level(q_position, parameters.side_length))
                upper = parameters.scale * complex(RRCFilters._qam_axis_level(i_position, parameters.side_length), RRCFilters._qam_axis_level(q_position + 1, parameters.side_length))
                @test sum(qamdemod([lower], M; gray=true) .!= qamdemod([upper], M; gray=true)) == 1
            end
        end

        for M in (8, 32, 128, 512), gray in (true, false)
            table = RRCFilters._qam_nonsquare_table(M, gray)
            @test length(table) == M
            @test length(unique(table)) == M
            max_level = maximum(level -> max(abs(level[1]), abs(level[2])), table)
            @test !any(table) do (i_level, q_level)
                abs(i_level) == max_level && abs(q_level) == max_level
            end
        end

        let rng = Xoshiro(999)
            for M in (8, 32, 128, 512)
                parameters = RRCFilters._qam_parameters(M)
                bin_table = RRCFilters._qam_nonsquare_table(M, false)
                for _ in 1:20_000
                    i_value = (rand(rng) - 0.5) * 2 * (parameters.i_side + 4)
                    q_value = (rand(rng) - 0.5) * 2 * (parameters.q_side + 4)
                    fast_label = RRCFilters._qam_nonsquare_fast_bin_label(i_value, q_value, parameters)
                    isnothing(fast_label) && continue
                    @test fast_label == RRCFilters._qam_nonsquare_nearest_label(i_value, q_value, bin_table)
                end
            end
        end

        for (M, snr_db, seed) in ((4, 30.0, 31), (16, 35.0, 32), (64, 40.0, 33), (256, 45.0, 34), (8, 28.0, 35), (32, 36.0, 36), (128, 42.0, 37), (512, 48.0, 38))
            bits = repeat(qam_label_bits(M), 16)
            symbols = qammod(bits, M)
            @test mean(abs2, symbols) ≈ 1.0 rtol=64eps(Float64)
            @test qamdemod(awgn(Xoshiro(seed), symbols, snr_db), M) == bits
        end

        @test isempty(qammod(Int[], 4))
        @test isempty(qamdemod(ComplexF64[], 4))
        @test isempty(qammod(Int[], 32))
        @test isempty(qamdemod(ComplexF64[], 32))
        @test qamdemod(ComplexF64[complex(floatmax(Float64), floatmax(Float64))], 32; gray=false) ==
              BitVector([true, true, false, false, false])
        zero_based_bits = ZeroBasedVector([false, false, true, true])
        @test qammod(zero_based_bits, 16) == qammod(zero_based_bits.data, 16)
        @test_throws ArgumentError qammod([0, 1, 0], 4)
        @test_throws ArgumentError qammod([0, 2], 4)
        @test_throws ArgumentError qammod([0, 1], 8)
        @test_throws ArgumentError qammod([0, 1], true)
        @test_throws ArgumentError qammod([0, 1], 2)
        @test_throws ArgumentError qammod([0, 1], 12)
        @test_throws ArgumentError qammod([0, 1], 2048)
        @test_throws ArgumentError qammod([0, 1], typemax(Int))
        @test_throws MethodError qammod([0, 1], 4.0)
        @test_throws ArgumentError qamdemod(ComplexF64[NaN + 0im], 4)
        @test_throws ArgumentError qamdemod(ComplexF64[Inf + 0im], 4)
        @test_throws ArgumentError qamdemod(ComplexF64[0 + Inf * im], 4)
        @test_throws MethodError qamdemod([1.0], 4)

        scale16 = RRCFilters._qam_parameters(16).scale
        @test qamdemod(ComplexF64[0 + 0im], 16; gray=true) == BitVector([true, true, false, true])
        @test qamdemod(ComplexF64[scale16 * (-2 + 3im)], 16; gray=false) == BitVector([false, true, false, false])
        @test qamdemod(ComplexF64[100 + 100im], 16; gray=false) == BitVector([true, true, false, false])
        @test qamdemod(ComplexF64[-100 - 100im], 16; gray=false) == BitVector([false, false, true, true])
        @test @inferred(qammod(BitVector([false, true, true, false]), 4)) isa Vector{ComplexF64}
        @test @inferred(qamdemod(ComplexF64[1 + im], 4)) isa BitVector
        @test isempty(Test.detect_ambiguities(RRCFilters))
    end

    @testset "PSK helper foundations" begin
        @test RRCFilters._psk_parameters(2) == (M=2, bits_per_symbol=1)
        @test RRCFilters._psk_parameters(16).bits_per_symbol == 4
        @test_throws ArgumentError RRCFilters._psk_parameters(true)
        @test_throws ArgumentError RRCFilters._psk_parameters(1)
        @test_throws ArgumentError RRCFilters._psk_parameters(0)
        @test_throws ArgumentError RRCFilters._psk_parameters(-2)
        @test_throws ArgumentError RRCFilters._psk_parameters(12)
        @test_throws ArgumentError RRCFilters._psk_parameters(typemax(Int))

        for M in (2, 4, 8, 16, 32), position in 0:(M - 1)
            @test RRCFilters._psk_nearest_position(cis(2pi * position / M), M, 0.0) == position
        end
    end

    @testset "PSK modulation and demodulation" begin
        function psk_label_bits(M)
            parameters = RRCFilters._psk_parameters(M)
            bits = falses(M * parameters.bits_per_symbol)
            for label in 0:(M - 1)
                RRCFilters._integer_to_bits!(bits, label * parameters.bits_per_symbol + 1, parameters.bits_per_symbol, label)
            end
            return bits
        end

        for M in (2, 4, 8, 16, 32), gray in (true, false)
            parameters = RRCFilters._psk_parameters(M)
            bits = psk_label_bits(M)
            original_bits = copy(bits)
            symbols = pskmod(bits, M; gray)

            @test length(symbols) == M
            @test length(unique(symbols)) == M
            @test all(symbol -> isapprox(abs(symbol), 1.0; atol=8eps(Float64)), symbols)
            @test mean(abs2, symbols) ≈ 1.0 rtol=8eps(Float64)
            @test pskdemod(symbols, M; gray) == bits
            @test pskmod(@view(bits[:]), M; gray) == symbols
            @test pskdemod(@view(symbols[:]), M; gray) == bits
            @test bits == original_bits

            reference_bits = bits[1:parameters.bits_per_symbol]
            reference_symbol = pskmod(reference_bits, M; gray)[1]
            perturbed = [reference_symbol * cis((2pi / M) / 4)]
            @test pskdemod(perturbed, M; gray) == reference_bits
        end

        for M in (2, 4, 8, 16, 32)
            for position in 0:(M - 1)
                left = ComplexF64[cis(2pi * position / M)]
                right = ComplexF64[cis(2pi * mod(position + 1, M) / M)]
                @test sum(pskdemod(left, M; gray=true) .!= pskdemod(right, M; gray=true)) == 1
            end
        end

        let rng = Xoshiro(1999)
            for M in (2, 4, 8, 16), phase_offset in (0.0, 0.3, -0.3, pi/4, pi/2, pi, -pi, 2pi, 100 * 2pi, -100 * 2pi, 12.345)
                rotation = cis(-phase_offset)
                rotation8 = cis(-phase_offset - pi / 8)
                rotation16 = cis(-phase_offset - pi / 16)
                for _ in 1:50_000
                    magnitude = 0.1 + 5.0 * rand(rng)
                    symbol = magnitude * cis(2pi * rand(rng))
                    fast = M == 8 ? RRCFilters._psk_octant_position(symbol, rotation8) :
                           M == 16 ? RRCFilters._psk_hexadecant_position(symbol, rotation16) :
                           RRCFilters._psk_fast_position(symbol, M, rotation)
                    isnothing(fast) && continue
                    @test fast == RRCFilters._psk_nearest_position(symbol, M, phase_offset)
                end
                for k in 0:(M - 1), perturb in (0.0, 1e-15, -1e-15, 1e-12, -1e-12, 1e-7, -1e-7)
                    boundary_angle = phase_offset + (pi / M) + k * (2pi / M) + perturb
                    symbol = cis(boundary_angle)
                    @test RRCFilters._psk_demodulate_position(symbol, M, phase_offset, rotation, rotation8, rotation16) == RRCFilters._psk_nearest_position(symbol, M, phase_offset)
                end
                for k in 0:(M - 1), perturb in (0.0, 1e-15, -1e-15, 1e-12, -1e-12, 1e-7, -1e-7)
                    center_angle = phase_offset + k * (2pi / M) + perturb
                    symbol = cis(center_angle)
                    @test RRCFilters._psk_demodulate_position(symbol, M, phase_offset, rotation, rotation8, rotation16) == RRCFilters._psk_nearest_position(symbol, M, phase_offset)
                end
                @test RRCFilters._psk_demodulate_position(ComplexF64(0, 0), M, phase_offset, rotation, rotation8, rotation16) == RRCFilters._psk_nearest_position(ComplexF64(0, 0), M, phase_offset)
            end
        end

        for (M, snr_db, seed) in ((2, 20.0, 41), (4, 25.0, 42), (8, 32.0, 43), (16, 38.0, 44))
            bits = repeat(psk_label_bits(M), 16)
            symbols = pskmod(bits, M)
            @test mean(abs2, symbols) ≈ 1.0 rtol=8eps(Float64)
            @test pskdemod(awgn(Xoshiro(seed), symbols, snr_db), M) == bits
        end

        bits8 = psk_label_bits(8)
        base_symbols = pskmod(bits8, 8; phase_offset=0.3)
        @test pskmod(bits8, 8; phase_offset=0.3 + 2pi) ≈ base_symbols rtol=1e-9
        @test pskmod(bits8, 8; phase_offset=0.3 + 100 * 2pi) ≈ base_symbols rtol=1e-9
        @test pskdemod(base_symbols, 8; phase_offset=0.3 + 2pi) == pskdemod(base_symbols, 8; phase_offset=0.3)
        @test pskdemod(base_symbols, 8; phase_offset=0.3 - 100 * 2pi) == pskdemod(base_symbols, 8; phase_offset=0.3)

        @test isempty(pskmod(Int[], 4))
        pskmod(Int[], 1 << 20)
        @test @allocated(pskmod(Int[], 1 << 20)) < 1024
        @test isempty(pskdemod(ComplexF64[], 4))
        zero_based_bits = ZeroBasedVector([false, false, true, true])
        @test pskmod(zero_based_bits, 4) == pskmod(zero_based_bits.data, 4)
        @test_throws ArgumentError pskmod([0, 1, 0], 4)
        @test_throws ArgumentError pskmod([0, 2], 4)
        @test_throws ArgumentError pskmod([0, 1], 3)
        @test_throws ArgumentError pskmod([0, 1], true)
        @test_throws ArgumentError pskmod([0, 1], 1)
        @test_throws ArgumentError pskmod([0, 1], 0)
        @test_throws ArgumentError pskmod([0, 1], -2)
        @test_throws ArgumentError pskmod([0, 1], 12)
        @test_throws ArgumentError pskmod([0, 1], typemax(Int))
        @test_throws MethodError pskmod([0, 1], 4.0)
        @test_throws ArgumentError pskmod([0, 1], 4; phase_offset=true)
        @test_throws ArgumentError pskmod([0, 1], 4; phase_offset=NaN)
        @test_throws ArgumentError pskmod([0, 1], 4; phase_offset=Inf)
        @test_throws ArgumentError pskdemod(ComplexF64[NaN + 0im], 4)
        @test_throws ArgumentError pskdemod(ComplexF64[Inf + 0im], 4)
        @test_throws ArgumentError pskdemod(ComplexF64[0 + Inf * im], 4)
        @test_throws ArgumentError pskdemod(ComplexF64[1 + 0im], 4; phase_offset=true)
        @test_throws ArgumentError pskdemod(ComplexF64[1 + 0im], 4; phase_offset=NaN)
        @test_throws MethodError pskdemod([1.0], 4)

        @test pskdemod(ComplexF64[0 + 0im], 4) == BitVector([false, false])
        @test @inferred(pskmod(BitVector([false, true]), 4)) isa Vector{ComplexF64}
        @test @inferred(pskdemod(ComplexF64[1 + 0im], 4)) isa BitVector
        @test isempty(Test.detect_ambiguities(RRCFilters))
    end

    @testset "upfirdn" begin
        @testset "input validation" begin
            @test_throws ArgumentError upfirdn([1.0], [1.0], true, 1)
            @test_throws ArgumentError upfirdn([1.0], [1.0], 1, true)
            @test_throws ArgumentError upfirdn([1.0], [1.0], 0, 1)
            @test_throws ArgumentError upfirdn([1.0], [1.0], 1, 0)
            @test_throws ArgumentError upfirdn([1.0], [1.0], -1, 1)
            @test_throws ArgumentError upfirdn([1.0], [1.0], 1, -1)
            @test_throws ArgumentError upfirdn([1.0], Float64[], 1, 1)
            @test_throws ArgumentError upfirdn([NaN], [1.0], 1, 1)
            @test_throws ArgumentError upfirdn([Inf], [1.0], 1, 1)
            @test_throws ArgumentError upfirdn([1.0], [NaN], 1, 1)
            @test_throws ArgumentError upfirdn([1.0], [Inf], 1, 1)
            @test_throws ArgumentError upfirdn(ComplexF64[1 + Inf * im], [1.0], 1, 1)
            @test_throws ArgumentError upfirdn([1.0], ComplexF64[1 + Inf * im], 1, 1)
            @test_throws MethodError upfirdn([1.0], [1.0], 1.0, 1)
            @test_throws MethodError upfirdn([1.0], [1.0], 1, 1.0)
        end

        @testset "independent BigFloat reference" begin
            rng = Xoshiro(2001)
            for _ in 1:200
                nx = rand(rng, 1:8)
                nh = rand(rng, 1:6)
                p = rand(rng, 1:4)
                q = rand(rng, 1:4)
                x = rand(rng, Bool) ? randn(rng, ComplexF64, nx) : randn(rng, nx)
                h = rand(rng, Bool) ? randn(rng, ComplexF64, nh) : randn(rng, nh)
                expected = reference_upfirdn(x, h, p, q)
                actual = upfirdn(x, h, p, q)
                @test length(actual) == length(expected)
                @test actual ≈ expected rtol=1e-9
            end
        end

        @testset "steady-region fast path, boundary math" begin
            # Exercises _upfirdn_steady_region/_upfirdn_boundary_range! across
            # nx sizes that land exactly on, just below, and well past the
            # steady-region threshold, for both shapes this library actually
            # uses (q_int==1 pure upsampling, p_int==1 pure decimation), for
            # both real and complex signals, and for a short (25-tap) and
            # long (513-tap) filter.
            rng = Xoshiro(2003)
            short_filter = rcosdesign(0.25, 6, 4; shape=:sqrt)  # 25 taps
            long_filter = rcosdesign(0.25, 32, 16; shape=:sqrt)  # 513 taps
            for filter in (short_filter, long_filter), p in (1, 2, 3, 4, 16, 32), q in (1, 3, 4, 16)
                p != 1 && q != 1 && continue  # general rational (both >1) isn't fast-pathed; skip here
                for nx in (1, 2, 3, 5, 8, 13, 21, 34, 55, 100, 1000)
                    x_real = randn(rng, nx)
                    x_complex = randn(rng, ComplexF64, nx)
                    for x in (x_real, x_complex)
                        expected = reference_upfirdn(x, filter, p, q)
                        actual = upfirdn(x, filter, p, q)
                        @test length(actual) == length(expected)
                        @test actual ≈ expected rtol=1e-9 atol=1e-12
                    end
                end
            end
        end

        @testset "degenerate and edge cases" begin
            x = [1.0, 2.0, 3.0]
            h = [1.0, 0.5]
            @test upfirdn(x, h) == [1.0, 2.5, 4.0, 1.5]
            @test upfirdn(x, h, 1, 1) == upfirdn(x, h)
            @test upfirdn(x, h) ≈ reference_upfirdn(x, h, 1, 1) rtol=1e-14
            @test isempty(upfirdn(Float64[], h, 3, 2))
            @test upfirdn(Float64[], h, 3, 2) isa Vector{Float64}
            @test isempty(upfirdn(ComplexF64[], h, 3, 2))
            @test upfirdn(ComplexF64[], h, 3, 2) isa Vector{ComplexF64}

            @test upfirdn([2.0], [1.0, -0.5], 1 << 20, 1) == [2.0, -1.0]
            upfirdn([1.0], [1.0], 1 << 20, 1)
            @test @allocated(upfirdn([1.0], [1.0], 1 << 20, 1)) < 1024
        end

        @testset "RRC self-convolution invariant" begin
            for beta in (0.0, 0.25, 0.5, 1.0), span in (2, 4, 6), sps in (2, 4)
                rrc = rcosdesign(beta, span, sps; shape=:sqrt)
                combined = upfirdn(rrc, rrc, 1, 1)
                @test length(combined) == 2 * span * sps + 1
                @test combined ≈ reverse(combined) rtol=1e-12 atol=1e-14
                @test argmax(combined) == span * sps + 1
            end
        end

        @testset "vectorized long-filter path" begin
            rng = Xoshiro(2002)
            long_filter = rcosdesign(0.25, 32, 16; shape=:sqrt)
            x = randn(rng, ComplexF64, 64)
            for (p, q) in ((16, 1), (1, 16), (3, 5))
                expected = reference_upfirdn(x, long_filter, p, q)
                actual = upfirdn(x, long_filter, p, q)
                @test actual ≈ expected rtol=1e-9
            end
        end

        @testset "types and immutability" begin
            real32 = Float32[1, -1, 0.5]
            real64 = Float64[1, -1, 0.5]
            complex32 = ComplexF32[1 + im, -1 - im]
            complex64 = ComplexF64[1 + im, -1 - im]
            bigfloats = BigFloat[big"1.0", big"-1.0"]
            integers = Int[1, -1, 2]
            filter32 = Float32[0.5, 0.25]
            filter64 = Float64[0.5, 0.25]

            @test eltype(upfirdn(real32, filter32, 2, 1)) === Float32
            @test eltype(upfirdn(real64, filter32, 2, 1)) === Float64
            @test eltype(upfirdn(real32, filter64, 2, 1)) === Float64
            @test eltype(upfirdn(complex32, filter32, 2, 1)) === ComplexF32
            @test eltype(upfirdn(complex64, filter32, 2, 1)) === ComplexF64
            @test eltype(upfirdn(real64, complex64, 2, 1)) === ComplexF64
            @test eltype(upfirdn(bigfloats, filter64, 2, 1)) === BigFloat
            @test eltype(upfirdn(integers, filter64, 2, 1)) === Float64

            original_x = copy(real64)
            original_h = copy(filter64)
            upfirdn(real64, filter64, 3, 2)
            @test real64 == original_x
            @test filter64 == original_h

            zero_based_x = ZeroBasedVector([1.0, -1.0, 2.0, 0.5])
            @test upfirdn(zero_based_x, filter64, 3, 1) == upfirdn(zero_based_x.data, filter64, 3, 1)
            @test upfirdn(@view(real64[:]), filter64, 2, 1) == upfirdn(real64, filter64, 2, 1)

            @test @inferred(upfirdn(real64, filter64, 2, 1)) isa Vector{Float64}
            @test isempty(Test.detect_ambiguities(RRCFilters))
        end
    end

    @testset "pulse shaping and matched filtering link" begin
        function shaped_link(bits, M, sps, span, beta, snr_db, rng)
            symbols = qammod(bits, M)
            txfilter = rcosdesign(beta, span, sps; shape=:sqrt)
            waveform = upfirdn(symbols, txfilter, sps)
            noisy_waveform = awgn(rng, waveform, snr_db)
            rxfilter = rcosdesign(beta, span, sps; shape=:sqrt)
            matched = upfirdn(noisy_waveform, rxfilter, 1, sps)
            symbol_count = length(symbols)
            return matched[(span + 1):(span + symbol_count)]
        end

        for (M, sps, span, beta, snr_db, seed) in ((4, 4, 6, 0.25, 25.0, 401), (16, 4, 6, 0.25, 32.0, 402))
            bits_per_symbol = trailing_zeros(M)
            bits = bitrand(Xoshiro(seed), 2048 * bits_per_symbol)
            ideal_symbols = qammod(bits, M)

            recovered_symbols = shaped_link(bits, M, sps, span, beta, snr_db, Xoshiro(seed + 1))
            @test length(recovered_symbols) == length(bits) ÷ bits_per_symbol
            recovered_bits = qamdemod(recovered_symbols, M)
            errors, ber = biterr(bits, recovered_bits)
            @test errors == 0
            @test ber == 0.0
            @test recovered_bits == bits

            # Isolates the residual ISI from RRC truncation (finite span) from
            # AWGN, at a per-sample SNR high enough to be effectively noiseless.
            near_noiseless_symbols = shaped_link(bits, M, sps, span, beta, 200.0, Xoshiro(seed + 2))
            @test evm(ideal_symbols, near_noiseless_symbols) < 5.0
        end
    end

    @testset "carriersync" begin
        @testset "input validation" begin
            seq = ComplexF64[1.0 + 0.0im, 0.0 + 1.0im]

            @test_throws ArgumentError carriersync(seq; modulation=:bogus)
            @test_throws ArgumentError carriersync(seq; sps=0)
            @test_throws ArgumentError carriersync(seq; sps=-1)
            @test_throws ArgumentError carriersync(seq; sps=true)
            @test_throws TypeError carriersync(seq; sps=1.5)
            @test_throws ArgumentError carriersync(seq; damping_factor=0.0)
            @test_throws ArgumentError carriersync(seq; damping_factor=-1.0)
            @test_throws ArgumentError carriersync(seq; damping_factor=true)
            @test_throws ArgumentError carriersync(seq; damping_factor=NaN)
            @test_throws ArgumentError carriersync(seq; damping_factor=Inf)
            @test_throws ArgumentError carriersync(seq; loop_bandwidth=0.0)
            @test_throws ArgumentError carriersync(seq; loop_bandwidth=-0.1)
            @test_throws ArgumentError carriersync(seq; loop_bandwidth=1.1)
            @test_throws ArgumentError carriersync(seq; loop_bandwidth=true)
            @test_throws ArgumentError carriersync(seq; loop_bandwidth=NaN)
            @test_throws ArgumentError carriersync(seq; loop_bandwidth=Inf)
            @test_throws ArgumentError carriersync(ComplexF64[NaN + 0im])
            @test_throws ArgumentError carriersync(ComplexF64[Inf + 0im])
            @test_throws ArgumentError carriersync(ComplexF64[0 + Inf * im])

            # Mode/config validation runs unconditionally, even for empty x
            # (mirrors upfirdn's _validate_upfirdn_arguments).
            @test_throws ArgumentError carriersync(ComplexF64[]; modulation=:bogus)
            @test_throws ArgumentError carriersync(ComplexF64[]; loop_bandwidth=0.0)
        end

        @testset "phase ambiguity and convergence" begin
            # An exact axis/diagonal-aligned point, repeated with no rotation
            # or noise, is already a fixed point of every PED formula (every
            # branch multiplies by sign(0)=0), so the loop must stay locked at
            # exactly zero phase forever -- an exact (not approximate) check.
            for (modulation, point) in ((:qam, (1.0 + 1.0im) / sqrt(2)), (:bpsk, 1.0 + 0.0im), (:psk8, 1.0 + 0.0im))
                x = fill(ComplexF64(point), 500)
                out, phest = carriersync(x; modulation, sps=1)
                @test all(iszero, phest)
                @test out == x
            end

            # A small constant per-sample rotation must be frequency-tracked:
            # the phase estimate's tail slope should match the true offset.
            for modulation in (:qam, :bpsk, :psk8)
                true_offset = 0.01
                base = modulation == :qam ? (1.0 + 1.0im) / sqrt(2) : 1.0 + 0.0im
                x = ComplexF64[base * cis(true_offset * (index - 1)) for index in 1:2000]
                _, phest = carriersync(x; modulation, sps=1)
                tail = 1801:2000
                slope = (phest[tail[end]] - phest[tail[1]]) / (tail[end] - tail[1])
                @test isapprox(slope, true_offset; atol=0.002)
            end
        end

        @testset "8PSK detector and noisy acquisition" begin
            rng = Xoshiro(771)
            for _ in 1:20_000
                sample = randn(rng) + im * randn(rng)
                actual = RRCFilters._carriersync_phase_error(sample, Val(:psk8))
                position = mod(round(Int, angle(sample) * 4 / pi, RoundNearestTiesAway), 8)
                expected = imag(sample * conj(cis(2pi * position / 8)))
                @test actual ≈ expected rtol=1e-12 atol=1e-12
            end
            for index in 0:7, offset in (-1e-6, 1e-6)
                sample = cis((2index + 1) * pi / 8 + offset)
                actual = RRCFilters._carriersync_phase_error(sample, Val(:psk8))
                position = mod(round(Int, angle(sample) * 4 / pi, RoundNearestTiesAway), 8)
                expected = imag(sample * conj(cis(2pi * position / 8)))
                @test actual ≈ expected rtol=1e-12 atol=1e-12
            end

            scenarios = ((1, 0.0005, 0.1), (2, 0.0015, -0.3), (3, -0.002, 0.7), (4, 0.003, -1.0),
                         (5, -0.001, 1.2), (6, 0.0025, -2.0), (7, 0.004, 0.4), (8, -0.003, -0.8))
            for (seed, frequency, phase) in scenarios
                bits = bitrand(Xoshiro(seed), 3 * 4096)
                symbols = pskmod(bits, 8)
                impaired = freqoffset(awgn(Xoshiro(seed + 100), symbols, 30.0);
                                      freq_offset=frequency / (2pi), phase_offset=phase)
                corrected, _ = carriersync(impaired; modulation=:psk8, sps=1)
                errors = minimum(0:7) do rotation
                    biterr(bits[901:end], pskdemod(corrected[301:end] .* cis(2pi * rotation / 8), 8))[1]
                end
                @test errors == 0
            end
        end

        @testset "types, immutability, and generic indexing" begin
            seq64 = ComplexF64[1.0 + 0.5im, -0.3 + 0.8im, 0.9 - 0.2im, -0.6 - 0.6im]
            seq32 = ComplexF32.(seq64)
            seqbig = Complex{BigFloat}.(seq64)
            original = copy(seq64)

            @test eltype(first(carriersync(seq64))) === ComplexF64
            @test eltype(last(carriersync(seq64))) === Float64
            @test eltype(first(carriersync(seq32))) === ComplexF32
            @test eltype(last(carriersync(seq32))) === Float32
            @test eltype(first(carriersync(seqbig))) === Complex{BigFloat}
            @test eltype(last(carriersync(seqbig))) === BigFloat
            @test eltype(first(carriersync(Int[1, -1, 0, 1]))) === ComplexF64
            @test eltype(first(carriersync(Rational{Int}[1 // 2, -1 // 2]))) === ComplexF64

            carriersync(seq64)
            @test seq64 == original

            zero_based = ZeroBasedVector(collect(seq64))
            @test carriersync(zero_based) == carriersync(zero_based.data)

            @test @inferred(carriersync(seq64)) isa Tuple{Vector{ComplexF64},Vector{Float64}}
        end

        @testset "degenerate and edge cases" begin
            out, phest = carriersync(ComplexF64[])
            @test isempty(out) && isempty(phest)
            @test out isa Vector{ComplexF64}
            @test phest isa Vector{Float64}

            single = ComplexF64[0.6 + 0.8im]
            out1, phest1 = carriersync(single)
            @test out1 == single
            @test phest1 == [0.0]
        end

        @testset "qammod/pskmod integration caveat" begin
            # Direct, deterministic confirmation of the geometric mismatch
            # itself, on _carriersync_phase_error, rather than relying on
            # emergent closed-loop behavior (which turned out to be more
            # forgiving in practice than the static geometry alone would
            # suggest -- see the closed-loop check below): a small rotation
            # error near a diagonal-aligned point gives a phase error
            # proportional to the rotation (a working restoring force); the
            # same small rotation near an axis-aligned point is dominated by
            # an O(1) term unrelated to the actual error, because exactly at
            # the axis (not perturbed) sign(0)=0 makes the error vanish only
            # in a non-generic, unstable way.
            epsilon = 0.01
            diagonal_error = RRCFilters._carriersync_phase_error(ComplexF64(cis(pi / 4 + epsilon)), Val(:qam))
            axis_error = RRCFilters._carriersync_phase_error(ComplexF64(cis(epsilon)), Val(:qam))
            @test abs(diagonal_error) < 0.1
            @test abs(axis_error) > 0.5

            # Closed-loop confirmation: qammod's own diagonal-aligned
            # constellation, and pskmod(bits, 4; phase_offset=pi/4) which
            # rotates onto the same diagonals, both correct perfectly under
            # :qam; pskmod(bits, 4) at its own default (axis-aligned, the
            # documented caveat) does not.
            bits = bitrand(Xoshiro(900), 2 * 2000)
            true_offset = 0.001

            function best_qam_errors(corrected, M, bits)
                return minimum(0:(M - 1)) do rotation
                    rotated = corrected .* cis(2pi * rotation / M)
                    errors, _ = biterr(bits, qamdemod(rotated, M))
                    return errors
                end
            end

            aligned_symbols = pskmod(bits, 4; phase_offset=pi / 4)
            misaligned_symbols = pskmod(bits, 4)
            aligned_offset = ComplexF64[s * cis(true_offset * (i - 1)) for (i, s) in enumerate(aligned_symbols)]
            misaligned_offset = ComplexF64[s * cis(true_offset * (i - 1)) for (i, s) in enumerate(misaligned_symbols)]

            aligned_corrected, _ = carriersync(aligned_offset; modulation=:qam, sps=1)
            misaligned_corrected, _ = carriersync(misaligned_offset; modulation=:qam, sps=1)

            @test best_qam_errors(aligned_corrected, 4, bits) == 0
            @test best_qam_errors(misaligned_corrected, 4, bits) > best_qam_errors(aligned_corrected, 4, bits)
        end
    end

    @testset "freqoffset" begin
        @testset "input validation" begin
            seq = ComplexF64[1.0 + 0.0im, 0.0 + 1.0im]

            @test_throws ArgumentError freqoffset(seq; freq_offset=true)
            @test_throws ArgumentError freqoffset(seq; freq_offset=NaN)
            @test_throws ArgumentError freqoffset(seq; freq_offset=Inf)
            @test_throws ArgumentError freqoffset(seq; phase_offset=true)
            @test_throws ArgumentError freqoffset(seq; phase_offset=NaN)
            @test_throws ArgumentError freqoffset(seq; phase_offset=-Inf)
            @test_throws ArgumentError freqoffset(seq; sample_rate=true)
            @test_throws ArgumentError freqoffset(seq; sample_rate=NaN)
            @test_throws ArgumentError freqoffset(seq; sample_rate=Inf)
            @test_throws ArgumentError freqoffset(seq; sample_rate=0.0)
            @test_throws ArgumentError freqoffset(seq; sample_rate=-1.0)
            @test_throws ArgumentError freqoffset(ComplexF64[NaN + 0im])
            @test_throws ArgumentError freqoffset(ComplexF64[Inf + 0im])
            @test_throws ArgumentError freqoffset(ComplexF64[0 + Inf * im])
            @test_throws ArgumentError freqoffset(Float64[NaN])

            @test_throws ArgumentError freqoffset(ComplexF64[]; sample_rate=0.0)
            @test_throws ArgumentError freqoffset(ComplexF64[]; freq_offset=NaN)
            @test_throws ArgumentError freqoffset(Float32[1.0]; freq_offset=floatmax(Float64))
        end

        @testset "closed-form oscillator reference" begin
            x = ComplexF64[1 + 0im, 0 + 1im, -0.5 + 0.25im, 0.75 - 0.2im]
            frequency = 17.25
            phase = -0.37
            sample_rate = 128.0
            expected = ComplexF64[
                sample * cis(2pi * frequency * (index - 1) / sample_rate + phase)
                for (index, sample) in enumerate(x)
            ]
            original = copy(x)
            @test freqoffset(x; freq_offset=frequency, phase_offset=phase, sample_rate) ≈ expected rtol=2e-14 atol=2e-14
            @test x == original

            for target in (0.0, 1 // 8, -3 // 7, 0.1, pi, 1e-12, 1e6 + 1 // 3)
                numerator, denominator = RRCFilters._freqoffset_rat(Float64(target), eps(Float64))
                @test isapprox(numerator / denominator, Float64(target); atol=eps(Float64) * max(1, abs(Float64(target))))
            end
        end

        @testset "carriersync round trip" begin
            bits = bitrand(Xoshiro(2028), 2 * 2000)
            symbols = pskmod(bits, 4; phase_offset=pi / 4)
            rad_per_sample = 0.01
            impaired = freqoffset(symbols; freq_offset=rad_per_sample / (2pi), phase_offset=0.4)
            _, phest = carriersync(impaired; modulation=:qam, sps=1)
            tail = 1801:2000
            slope = (phest[tail[end]] - phest[tail[1]]) / (tail[end] - tail[1])
            @test isapprox(slope, rad_per_sample; atol=0.002)
        end

        @testset "types, immutability, and generic indexing" begin
            seq64 = ComplexF64[1.0 + 0.5im, -0.3 + 0.8im, 0.9 - 0.2im, -0.6 - 0.6im]
            seq32 = ComplexF32.(seq64)
            seqbig = Complex{BigFloat}.(seq64)
            original = copy(seq64)
            kwargs = (freq_offset=0.01, phase_offset=0.3)

            @test eltype(freqoffset(seq64; kwargs...)) === ComplexF64
            @test eltype(freqoffset(seq32; kwargs...)) === ComplexF32
            @test eltype(freqoffset(seqbig; kwargs...)) === Complex{BigFloat}
            @test eltype(freqoffset(Float64[1.0, -0.5]; kwargs...)) === ComplexF64
            @test eltype(freqoffset(Float32[1.0, -0.5]; kwargs...)) === ComplexF32
            @test eltype(freqoffset(Int[1, -1, 0]; kwargs...)) === ComplexF64
            @test eltype(freqoffset(Rational{Int}[1 // 2, -1 // 2]; kwargs...)) === ComplexF64

            freqoffset(seq64; kwargs...)
            @test seq64 == original

            zero_based = ZeroBasedVector(collect(seq64))
            @test freqoffset(zero_based; kwargs...) == freqoffset(zero_based.data; kwargs...)

            @test @inferred(freqoffset(seq64; kwargs...)) isa Vector{ComplexF64}
        end

        @testset "degenerate and edge cases" begin
            @test freqoffset(ComplexF64[]) isa Vector{ComplexF64}
            @test isempty(freqoffset(ComplexF64[]))
            @test freqoffset(Float64[]) isa Vector{ComplexF64}
            @test freqoffset(Float32[]) isa Vector{ComplexF32}
            @test freqoffset(Complex{BigFloat}[]) isa Vector{Complex{BigFloat}}

            seq = ComplexF64[1.0 + 0.5im, -0.3 + 0.8im, 0.0 + 0.0im]
            @test freqoffset(seq) == seq
            @test freqoffset(Float64[1.0, -2.5, 0.0]) == ComplexF64[1.0, -2.5, 0.0]

            @test freqoffset(Float64[1.0]; phase_offset=pi / 2) == ComplexF64[im]

            @test freqoffset(Float64[1.0, 1.0, 1.0, 1.0]; freq_offset=0.5) ==
                  ComplexF64[1.0, -1.0, 1.0, -1.0]
            @test freqoffset(Float64[1.0]; freq_offset=floatmax(Float64), sample_rate=floatmin(Float64)) == ComplexF64[1.0]
        end
    end

    @testset "symbolsync" begin
        @testset "input validation" begin
            seq = ComplexF64[1.0 + 0.0im, 0.0 + 1.0im, -1.0 + 0.0im, 0.0 - 1.0im]

            @test_throws UndefKeywordError symbolsync(seq)
            @test_throws ArgumentError symbolsync(seq; sps=4, timing_error_detector=:bogus)
            @test_throws ArgumentError symbolsync(seq; sps=0)
            @test_throws ArgumentError symbolsync(seq; sps=1)
            @test_throws ArgumentError symbolsync(seq; sps=-1)
            @test_throws ArgumentError symbolsync(seq; sps=true)
            @test_throws TypeError symbolsync(seq; sps=1.5)
            @test_throws ArgumentError symbolsync(seq; sps=4, damping_factor=0.0)
            @test_throws ArgumentError symbolsync(seq; sps=4, damping_factor=-1.0)
            @test_throws ArgumentError symbolsync(seq; sps=4, damping_factor=true)
            @test_throws ArgumentError symbolsync(seq; sps=4, damping_factor=NaN)
            @test_throws ArgumentError symbolsync(seq; sps=4, damping_factor=Inf)
            @test_throws ArgumentError symbolsync(seq; sps=4, loop_bandwidth=0.0)
            @test_throws ArgumentError symbolsync(seq; sps=4, loop_bandwidth=-0.1)
            @test_throws ArgumentError symbolsync(seq; sps=4, loop_bandwidth=1.0)
            @test_throws ArgumentError symbolsync(seq; sps=4, loop_bandwidth=1.1)
            @test_throws ArgumentError symbolsync(seq; sps=4, loop_bandwidth=true)
            @test_throws ArgumentError symbolsync(seq; sps=4, loop_bandwidth=NaN)
            @test_throws ArgumentError symbolsync(seq; sps=4, loop_bandwidth=Inf)
            @test_throws ArgumentError symbolsync(seq; sps=4, detector_gain=0.0)
            @test_throws ArgumentError symbolsync(seq; sps=4, detector_gain=-1.0)
            @test_throws ArgumentError symbolsync(seq; sps=4, detector_gain=true)
            @test_throws ArgumentError symbolsync(seq; sps=4, detector_gain=NaN)
            @test_throws ArgumentError symbolsync(seq; sps=4, detector_gain=Inf)
            @test_throws ArgumentError symbolsync(ComplexF64[NaN + 0im]; sps=4)
            @test_throws ArgumentError symbolsync(ComplexF64[Inf + 0im]; sps=4)
            @test_throws ArgumentError symbolsync(ComplexF64[0 + Inf * im]; sps=4)

            @test_throws ArgumentError symbolsync(ComplexF64[]; sps=4, timing_error_detector=:bogus)
            @test_throws ArgumentError symbolsync(ComplexF64[]; sps=4, loop_bandwidth=0.0)
        end

        @testset "TED and trigger unit tests" begin
            for ted in (:zero_crossing, :gardner, :early_late, :mueller_muller)
                tag = Val(ted)
                buf = Complex{Float64}[1.0, 0.5, -0.5, 1.0]
                mids = (buf[2] + buf[3]) / 2
                @test RRCFilters._symbolsync_ted_error(buf, mids, ComplexF64(0.2), tag) isa Float64

                history_one = fill(false, 4)
                history_one[1] = true
                @test RRCFilters._symbolsync_trigger(true, history_one, tag) == (ted != :early_late)

                history_prev = fill(false, 4)
                history_prev[end] = true
                @test RRCFilters._symbolsync_trigger(false, history_prev, tag) == (ted == :early_late)
            end
        end

        @testset "types, immutability, and generic indexing" begin
            seq64 = ComplexF64[1.0 + 0.5im, -0.3 + 0.8im, 0.9 - 0.2im, -0.6 - 0.6im, 0.2 + 0.1im, -0.4 - 0.3im]
            seq32 = ComplexF32.(seq64)
            seqbig = Complex{BigFloat}.(seq64)
            original = copy(seq64)

            @test eltype(first(symbolsync(seq64; sps=2))) === ComplexF64
            @test eltype(last(symbolsync(seq64; sps=2))) === Float64
            @test eltype(first(symbolsync(seq32; sps=2))) === ComplexF32
            @test eltype(last(symbolsync(seq32; sps=2))) === Float32
            @test eltype(first(symbolsync(seqbig; sps=2))) === Complex{BigFloat}
            @test eltype(last(symbolsync(seqbig; sps=2))) === BigFloat
            @test eltype(first(symbolsync(Int[1, -1, 0, 1, -1, 1]; sps=2))) === ComplexF64
            @test eltype(first(symbolsync(Rational{Int}[1 // 2, -1 // 2, 1 // 2, -1 // 2, 1 // 2, -1 // 2]; sps=2))) === ComplexF64

            symbolsync(seq64; sps=2)
            @test seq64 == original

            zero_based = ZeroBasedVector(collect(seq64))
            @test symbolsync(zero_based; sps=2) == symbolsync(zero_based.data; sps=2)

            @test @inferred(symbolsync(seq64; sps=2)) isa Tuple{Vector{ComplexF64},Vector{Float64}}
        end

        @testset "degenerate and edge cases" begin
            symbols, timing_error = symbolsync(ComplexF64[]; sps=4)
            @test isempty(symbols) && isempty(timing_error)
            @test symbols isa Vector{ComplexF64}
            @test timing_error isa Vector{Float64}

            single = ComplexF64[0.6 + 0.8im]
            symbols1, timing_error1 = symbolsync(single; sps=4)
            @test isempty(symbols1)
            @test timing_error1 == [0.0]
        end
    end

    @testset "timingoffset" begin
        @testset "input validation" begin
            seq = Float64[1.0, 2.0, 3.0, 4.0]

            @test_throws ArgumentError timingoffset(seq; delay=true)
            @test_throws ArgumentError timingoffset(seq; delay=NaN)
            @test_throws ArgumentError timingoffset(seq; delay=Inf)
            @test_throws ArgumentError timingoffset(seq; delay=-Inf)
        end

        @testset "integer-delay exactness" begin
            x = Float64.(1:10)
            for d in (-3, -1, 0, 1, 2, 5)
                y = timingoffset(x; delay=d)
                expected = zeros(Float64, length(x))
                for i in eachindex(x)
                    if i > 2
                        src = i - d - 2
                        1 <= src <= length(x) && (expected[i] = x[src])
                    end
                end
                @test y == expected
            end
        end

        @testset "fractional-delay phase shift on a sinusoid" begin
            n = 2000
            w = 0.03
            x = ComplexF64[cis(w * k) for k in 0:(n - 1)]
            for d in (0.0, 0.25, 0.5, 0.9, -1.3)
                y = timingoffset(x; delay=d)
                interior = 50:(n - 50)
                expected = ComplexF64[cis(w * (k - 1 - d - 2)) for k in interior]
                @test maximum(abs.(y[interior] .- expected)) < 0.08
            end
        end

        @testset "boundary zero-padding and energy sanity" begin
            x = Float64.(1:100)
            y = timingoffset(x; delay=3.4)
            @test all(iszero, y[1:3])
            @test !iszero(y[4])
            interior_energy_x = sum(abs2, x[10:90])
            interior_energy_y = sum(abs2, y[15:95])
            @test isapprox(interior_energy_y, interior_energy_x; rtol=0.05)
        end

        @testset "symbolsync round trip" begin
            bits = bitrand(Xoshiro(3031), 2 * 4000)
            symbols = pskmod(bits, 4; phase_offset=pi / 4)
            sps = 4
            span = 8
            txfilter = rcosdesign(0.25, span, sps; shape=:sqrt)
            shaped = upfirdn(symbols, txfilter, sps)
            true_delay = 1.7
            delayed = timingoffset(shaped; delay=true_delay)
            noisy = awgn(Xoshiro(3032), delayed, 30.0)
            rxfilter = rcosdesign(0.25, span, sps; shape=:sqrt)
            matched = upfirdn(noisy, rxfilter, 1, 1)

            recovered_symbols, timing_error = symbolsync(matched; sps=sps, timing_error_detector=:gardner)
            tail = timing_error[(end - 399):end]
            @test maximum(tail) - minimum(tail) < 0.05

            recovered_tail = recovered_symbols[(end - 499):end]
            demodulated = pskdemod(recovered_tail, 4; phase_offset=pi / 4)
            remodulated = pskmod(demodulated, 4; phase_offset=pi / 4)
            @test evm(remodulated, recovered_tail; mode=:percent) < 20.0
        end

        @testset "types, immutability, and generic indexing" begin
            seq64 = Float64[1.0, -0.5, 0.9, -0.6, 0.2, -0.4]
            seq32 = Float32.(seq64)
            seqbig = BigFloat.(seq64)
            cseq64 = ComplexF64.(seq64)
            original = copy(seq64)

            @test eltype(timingoffset(seq64; delay=0.3)) === Float64
            @test eltype(timingoffset(seq32; delay=0.3)) === Float32
            @test eltype(timingoffset(seqbig; delay=0.3)) === BigFloat
            @test eltype(timingoffset(cseq64; delay=0.3)) === ComplexF64
            @test eltype(timingoffset(Int[1, -1, 0, 1, -1, 1]; delay=0.3)) === Float64
            @test eltype(timingoffset(Rational{Int}[1 // 2, -1 // 2, 1 // 2, -1 // 2, 1 // 2, -1 // 2]; delay=0.3)) === Float64

            timingoffset(seq64; delay=0.3)
            @test seq64 == original

            zero_based = ZeroBasedVector(collect(seq64))
            @test timingoffset(zero_based; delay=0.3) == timingoffset(zero_based.data; delay=0.3)

            @test @inferred(timingoffset(seq64; delay=0.3)) isa Vector{Float64}
        end

        @testset "degenerate and edge cases" begin
            @test timingoffset(Float64[]) isa Vector{Float64}
            @test isempty(timingoffset(Float64[]))
            @test timingoffset(ComplexF64[]) isa Vector{ComplexF64}
            @test timingoffset(Float32[]) isa Vector{Float32}

            @test timingoffset(Float64[1.0, 2.0, 3.0]) == Float64[0.0, 0.0, 1.0]
        end
    end

    @testset "agc" begin
        @testset "input validation" begin
            seq = ComplexF64[1.0 + 0im, 2.0 + 0im]

            @test_throws ArgumentError agc(seq; desired_amplitude=true)
            @test_throws ArgumentError agc(seq; desired_amplitude=0.0)
            @test_throws ArgumentError agc(seq; desired_amplitude=-1.0)
            @test_throws ArgumentError agc(seq; desired_amplitude=NaN)
            @test_throws ArgumentError agc(seq; desired_amplitude=Inf)

            @test_throws ArgumentError agc(seq; step_size=true)
            @test_throws ArgumentError agc(seq; step_size=0.0)
            @test_throws ArgumentError agc(seq; step_size=-0.1)
            @test_throws ArgumentError agc(seq; step_size=NaN)
            @test_throws ArgumentError agc(seq; step_size=Inf)
            @test_throws ArgumentError agc(seq; step_size=1.0000001)
            @test_throws ArgumentError agc(seq; step_size=2.0)
            @test agc(seq; step_size=1.0)[1] isa Vector{ComplexF64}

            @test_throws ArgumentError agc(seq; initial_gain=true)
            @test_throws ArgumentError agc(seq; initial_gain=0.0)
            @test_throws ArgumentError agc(seq; initial_gain=-1.0)
            @test_throws ArgumentError agc(seq; initial_gain=NaN)
            @test_throws ArgumentError agc(seq; initial_gain=Inf)

            @test_throws ArgumentError agc(seq; min_gain=true)
            @test_throws ArgumentError agc(seq; min_gain=0.0)
            @test_throws ArgumentError agc(seq; min_gain=-1.0)
            @test_throws ArgumentError agc(seq; min_gain=NaN)
            @test_throws ArgumentError agc(seq; min_gain=Inf)

            @test_throws ArgumentError agc(seq; max_gain=true)
            @test_throws ArgumentError agc(seq; max_gain=0.0)
            @test_throws ArgumentError agc(seq; max_gain=-1.0)
            @test_throws ArgumentError agc(seq; max_gain=NaN)
            @test_throws ArgumentError agc(seq; max_gain=Inf)

            @test_throws ArgumentError agc(seq; min_gain=10.0, max_gain=1.0)
            @test_throws ArgumentError agc(seq; initial_gain=0.5, min_gain=1.0, max_gain=10.0)
            @test_throws ArgumentError agc(seq; initial_gain=100.0, min_gain=1.0, max_gain=10.0)

            @test_throws ArgumentError agc(ComplexF64[NaN + 0im, 1.0 + 0im])
            @test_throws ArgumentError agc(ComplexF64[Inf + 0im, 1.0 + 0im])
            @test_throws ArgumentError agc(ComplexF64[1.0 + NaN * im, 1.0 + 0im])
            @test_throws ArgumentError agc(ComplexF64[1.0 + Inf * im, 1.0 + 0im])
            @test_throws ArgumentError agc(Float64[NaN, 1.0])
            @test_throws ArgumentError agc(Float64[Inf, 1.0])

            @test_throws ArgumentError agc(Float64[]; step_size=0.0)
        end

        @testset "steady-state convergence and scale invariance" begin
            target_ratio = 1000.0
            convergence_indices = Int[]
            for amplitude in (1e-6, 1e-3, 1.0, 1e3, 1e6)
                desired = amplitude * target_ratio
                x = fill(amplitude, 2000)
                y, gain = agc(x; desired_amplitude=desired)
                @test gain[end] * amplitude ≈ desired rtol=1e-6
                idx = findfirst(i -> abs(gain[i] * amplitude / desired - 1.0) < 0.01, eachindex(gain))
                @test idx !== nothing
                push!(convergence_indices, idx)

                xc = fill(Complex(amplitude, 0.0), 2000)
                yc, gainc = agc(xc; desired_amplitude=desired)
                @test gainc[end] * amplitude ≈ desired rtol=1e-6
            end
            @test all(==(first(convergence_indices)), convergence_indices)
        end

        @testset "silence handling" begin
            x = zeros(500)
            y, gain = agc(x)
            @test all(iszero, y)
            @test all(isfinite, gain)
            @test gain[end] ≈ 1e4

            xc = zeros(ComplexF64, 500)
            yc, gainc = agc(xc)
            @test all(iszero, yc)
            @test all(isfinite, gainc)
            @test gainc[end] ≈ 1e4
        end

        @testset "extreme-amplitude and overflow-guard regression" begin
            x = fill(prevfloat(floatmax(Float64)), 50)
            y, gain = agc(x)
            @test all(isfinite, y)
            @test all(isfinite, gain)
            @test gain[end] ≈ 1e-4

            x1 = vcat([floatmax(Float64)], fill(1.0, 2000))
            y1, gain1 = agc(x1; initial_gain=1e4, max_gain=1e4)
            @test isfinite(y1[1])
            @test y1[1] == gain1[1] * x1[1]
            @test gain1[1] < 1e4
            @test isapprox(gain1[end], 1.0; atol=0.01)
        end

        @testset "types, immutability, and generic indexing" begin
            seq64 = Float64[1.0, -0.5, 0.9, -0.6, 0.2, -0.4, 0.7, -0.3]
            seq32 = Float32.(seq64)
            seq16 = Float16.(seq64)
            seqbig = BigFloat.(seq64)
            cseq64 = ComplexF64.(seq64)
            original = copy(seq64)

            @test eltype(agc(seq64)[1]) === Float64
            @test eltype(agc(seq64)[2]) === Float64
            @test eltype(agc(seq32)[1]) === Float32
            @test eltype(agc(seq32)[2]) === Float32
            @test eltype(agc(seq16)[1]) === Float16
            @test eltype(agc(seqbig)[1]) === BigFloat
            @test eltype(agc(cseq64)[1]) === ComplexF64
            @test eltype(agc(cseq64)[2]) === Float64
            @test eltype(agc(Int[1, -1, 2, -2, 3, -3, 4, -4])[1]) === Float64
            @test eltype(agc(Rational{Int}[1 // 2, -1 // 2, 1 // 2, -1 // 2, 1 // 2, -1 // 2, 1 // 2, -1 // 2])[1]) === Float64

            agc(seq64)
            @test seq64 == original

            zero_based = ZeroBasedVector(collect(seq64))
            @test agc(zero_based) == agc(zero_based.data)

            @test @inferred(agc(seq64)) isa Tuple{Vector{Float64},Vector{Float64}}
        end

        @testset "degenerate and edge cases" begin
            y, gain = agc(Float64[])
            @test isempty(y) && isempty(gain)
            @test eltype(y) === Float64 && eltype(gain) === Float64

            yc, gainc = agc(ComplexF64[])
            @test isempty(yc) && isempty(gainc)
            @test eltype(yc) === ComplexF64 && eltype(gainc) === Float64

            y32, gain32 = agc(Float32[])
            @test eltype(y32) === Float32 && eltype(gain32) === Float32

            y1, gain1 = agc([3.0])
            @test y1 == [3.0]
            @test gain1 == [1.0]

            y1b, gain1b = agc([3.0]; initial_gain=2.0, min_gain=2.0, max_gain=2.0)
            @test y1b == [6.0]
            @test gain1b == [2.0]

            xpin = [1.0, 5.0, 0.1, 10.0]
            ypin, gainpin = agc(xpin; step_size=0.5, initial_gain=2.0, min_gain=2.0, max_gain=2.0)
            @test all(==(2.0), gainpin)
            @test ypin == 2.0 .* xpin
        end

        @testset "QAM amplitude-swing round trip" begin
            M = 16
            bits_per_symbol = trailing_zeros(M)
            half = 6000
            symbol_count = 2 * half
            bits = bitrand(Xoshiro(9001), symbol_count * bits_per_symbol)
            symbols = qammod(bits, M)

            noisy = awgn(Xoshiro(9002), symbols, 25.0; signal_power=1.0)

            swung = copy(noisy)
            swung[1:half] .*= 0.05
            swung[(half + 1):end] .*= 20.0

            direct_demod = qamdemod(swung, M)
            _, direct_ber = biterr(bits, direct_demod)

            corrected, gain = agc(swung; desired_amplitude=1.0)
            corrected_demod = qamdemod(corrected, M)
            _, corrected_ber = biterr(bits, corrected_demod)

            @test direct_ber > 0.15
            @test corrected_ber < 0.02
            @test corrected_ber < direct_ber
        end
    end

    @testset "FEC" begin
        @testset "input validation" begin
            @test_throws ArgumentError poly2trellis(true, [7, 5])
            @test_throws ArgumentError poly2trellis(1, [7, 5])
            @test_throws ArgumentError poly2trellis(25, [7, 5])
            @test_throws ArgumentError poly2trellis(3, Int[])
            @test_throws ArgumentError poly2trellis(3, Bool[true, false])
            @test_throws ArgumentError poly2trellis(3, [-1, 5])
            @test_throws ArgumentError poly2trellis(3, [8, 5])
            @test_throws ArgumentError poly2trellis(3, fill(1, 17))
            t0 = poly2trellis(3, [0, 5])
            @test all(o -> o in (0, 1), t0.outputs)

            t = poly2trellis(3, [7, 5])
            @test_throws ArgumentError convenc([0, 1], t; initial_state=true)
            @test_throws ArgumentError convenc([0, 1], t; initial_state=-1)
            @test_throws ArgumentError convenc([0, 1], t; initial_state=4)
            @test_throws ArgumentError convenc([0, 2], t)
            @test convenc(Int[], t) == (BitVector(), 0)

            @test_throws ArgumentError vitdec([0, 1, 0], t, 1)
            @test isempty(vitdec(Int[], t, 1))
            @test_throws ArgumentError vitdec([0, 1, 0, 1], t, 1; mode=:cont)
            @test_throws ArgumentError vitdec([0, 1, 0, 1], t, true)
            @test_throws ArgumentError vitdec([0, 1, 0, 1], t, 0)
            @test_throws ArgumentError vitdec([0, 1, 0, 1], t, 3)
            @test_throws ArgumentError vitdec([0, 1, 0, 2], t, 1)

            @test_throws ArgumentError vitdec([0, 1, 0, 1], t, 1; decision_type=:foo)
            @test_throws ArgumentError vitdec([0, 1, 0, 1], t, 1; decision_type=:soft, num_soft_bits=true)
            @test_throws ArgumentError vitdec([0, 1, 0, 1], t, 1; decision_type=:soft, num_soft_bits=0)
            @test_throws ArgumentError vitdec([0, 1, 0, 1], t, 1; decision_type=:soft, num_soft_bits=-1)

            @test_throws ArgumentError vitdec([0, 1, -1, 1], t, 1; decision_type=:soft, num_soft_bits=3)
            @test_throws ArgumentError vitdec([0, 1, 8, 1], t, 1; decision_type=:soft, num_soft_bits=3)
            @test_throws ArgumentError vitdec([0.0, 1.0, 0.5, 1.0], t, 1; decision_type=:soft, num_soft_bits=3)
            @test (vitdec([0, 1, 7, 1], t, 1; decision_type=:soft, num_soft_bits=3); true)

            @test_throws ArgumentError vitdec([0.0, 1.0, NaN, -1.0], t, 1; decision_type=:unquant)
            @test_throws ArgumentError vitdec([0.0, 1.0, Inf, -1.0], t, 1; decision_type=:unquant)
            @test_throws ArgumentError vitdec(BigFloat[big"1e400", 1, -1, 1], t, 1; decision_type=:unquant)
            @test (vitdec([1.0, -1.0, 0.3, -0.8], t, 1; decision_type=:unquant); true)
            @test_throws ArgumentError vitdec([0, 1, 0, 1], t, 1; decision_type=:soft, num_soft_bits=53)

            invalid_states = copy(t.next_states)
            invalid_states[1, 1] = t.num_states
            @test_throws ArgumentError vitdec([0, 1], merge(t, (next_states=invalid_states,)), 1)

            invalid_outputs = copy(t.outputs)
            invalid_outputs[1, 1] = t.num_output_symbols
            @test_throws ArgumentError vitdec([0, 1], merge(t, (outputs=invalid_outputs,)), 1)
        end

        @testset "known trellis and encoder vector" begin
            trellis = poly2trellis(3, [7, 5])
            @test trellis.num_states == 4
            @test trellis.next_states == [0 2; 0 2; 1 3; 1 3]
            @test trellis.outputs == [0 3; 3 0; 2 1; 1 2]
            code, final_state = convenc(Bool[1, 0, 1, 1], trellis)
            @test code == Bool[1, 1, 1, 0, 0, 0, 0, 1]
            @test final_state == 3
        end

        @testset "num_soft_bits=1 exactly matches :hard" begin
            rng = Xoshiro(5151)
            trellis = poly2trellis(3, [7, 5])
            mem = trellis.constraint_length - 1
            msg = bitrand(rng, 100)
            msg_term = vcat(msg, falses(mem))
            code_term, _ = convenc(msg_term, trellis)
            noisy = copy(code_term)
            noisy[1:5:end] .= .!noisy[1:5:end]
            hard_decoded = vitdec(noisy, trellis, 20; mode=:term)
            soft_decoded = vitdec(Int.(noisy), trellis, 20; mode=:term, decision_type=:soft, num_soft_bits=1)
            @test soft_decoded == hard_decoded
        end

        @testset "unquant decoding realizes real coding gain over hard-decision" begin
            trellis = poly2trellis(7, [0o171, 0o133])
            mem = trellis.constraint_length - 1
            rng = Xoshiro(9090)
            msg = bitrand(rng, 4_000)
            msg_term = vcat(msg, falses(mem))
            coded, _ = convenc(msg_term, trellis)

            bpsk = ComplexF64[b ? -1.0 : 1.0 for b in coded]
            noisy = awgn(rng, bpsk, 1.0)
            real_samples = real.(noisy)

            hard_bits = real_samples .< 0
            hard_decoded = vitdec(hard_bits, trellis, 30; mode=:term)
            _, hard_ber = biterr(msg_term, hard_decoded)

            unquant_decoded = vitdec(real_samples, trellis, 30; mode=:term, decision_type=:unquant)
            _, unquant_ber = biterr(msg_term, unquant_decoded)

            @test unquant_ber < hard_ber
        end

        @testset "unquant decoding is scale-safe" begin
            trellis = poly2trellis(3, [7, 5])
            message = Bool[1, 0, 1, 1]
            code, _ = convenc(message, trellis)
            samples = Float64[bit ? -floatmax(Float64) : floatmax(Float64) for bit in code]
            @test vitdec(samples, trellis, 1; decision_type=:unquant) == message
        end

        @testset "round trips" begin
            trellis = poly2trellis(3, [7, 5])
            mem = trellis.constraint_length - 1
            rng = Xoshiro(4242)
            msg = bitrand(rng, 200)

            msg_term = vcat(msg, falses(mem))
            code_term, _ = convenc(msg_term, trellis)
            decoded_term = vitdec(code_term, trellis, 20; mode=:term)
            @test decoded_term == msg_term

            code_trunc, _ = convenc(msg, trellis)
            decoded_trunc = vitdec(code_trunc, trellis, 20; mode=:trunc)
            @test decoded_trunc == msg

            noisy_within = copy(code_term)
            noisy_within[15] = !noisy_within[15]
            @test vitdec(noisy_within, trellis, 20; mode=:term) == msg_term

            noisy_beyond = copy(code_term)
            noisy_beyond[1:3:end] .= .!noisy_beyond[1:3:end]
            decoded_beyond = vitdec(noisy_beyond, trellis, 20; mode=:term)
            @test length(decoded_beyond) == length(msg_term)
        end

        @testset "K=9 large-state round trip" begin
            trellis = poly2trellis(9, [0o561, 0o753])
            mem = trellis.constraint_length - 1
            rng = Xoshiro(777)
            msg = bitrand(rng, 500)
            msg_term = vcat(msg, falses(mem))
            code, _ = convenc(msg_term, trellis)
            decoded = vitdec(code, trellis, 30; mode=:term)
            @test decoded == msg_term
        end

        @testset "types, immutability, and generic indexing" begin
            trellis = poly2trellis(3, [7, 5])
            bits = Int[0, 1, 1, 0, 0, 1, 0, 1]
            original_bits = copy(bits)

            code, final_state = @inferred convenc(bits, trellis)
            @test bits == original_bits
            @test code isa BitVector

            original_code = copy(code)
            decoded = @inferred vitdec(code, trellis, 4)
            @test code == original_code
            @test decoded isa BitVector

            zero_bits = ZeroBasedVector(collect(bits))
            code_zb, final_state_zb = convenc(zero_bits, trellis)
            @test code_zb == code
            @test final_state_zb == final_state

            zero_code = ZeroBasedVector(collect(code))
            @test vitdec(zero_code, trellis, 4) == decoded

            soft_code = Int.(code)
            original_soft_code = copy(soft_code)
            decoded_soft = @inferred vitdec(soft_code, trellis, 4; decision_type=:soft, num_soft_bits=1)
            @test soft_code == original_soft_code
            @test decoded_soft isa BitVector
            @test decoded_soft == decoded

            unquant_code = Float64[b ? -1.0 : 1.0 for b in code]
            original_unquant_code = copy(unquant_code)
            decoded_unquant = @inferred vitdec(unquant_code, trellis, 4; decision_type=:unquant)
            @test unquant_code == original_unquant_code
            @test decoded_unquant isa BitVector
            @test decoded_unquant == decoded

            zero_soft_code = ZeroBasedVector(collect(soft_code))
            @test vitdec(zero_soft_code, trellis, 4; decision_type=:soft, num_soft_bits=1) == decoded_soft

            zero_unquant_code = ZeroBasedVector(collect(unquant_code))
            @test vitdec(zero_unquant_code, trellis, 4; decision_type=:unquant) == decoded_unquant
        end

        @testset "degenerate and edge cases" begin
            trellis = poly2trellis(3, [7, 5])
            @test convenc(Int[], trellis) == (BitVector(), 0)
            @test isempty(vitdec(Int[], trellis, 1))

            single_code, single_state = convenc([1], trellis)
            @test length(single_code) == 2
            @test vitdec(single_code, trellis, 1; mode=:trunc) == [true]
        end

        @test isempty(Test.detect_ambiguities(RRCFilters))
    end

    @testset "CRC" begin
        crc_round_trip_ok(msg, cfg) = begin
            codeword = crcgenerate(msg, cfg)
            decoded_msg, err = crcdetect(codeword, cfg)
            decoded_msg == BitVector(msg) && !any(err)
        end

        @testset "crcconfig validation" begin
            @test_throws ArgumentError crcconfig(Int[])
            @test_throws ArgumentError crcconfig([16, 12, 5, 5])
            @test_throws ArgumentError crcconfig([5, 12, 16, 0])
            @test_throws ArgumentError crcconfig([-1, 0])
            @test_throws ArgumentError crcconfig(Bool[true, false, false, true])
            @test_throws ArgumentError crcconfig([0, 1, 0, 1])
            @test_throws ArgumentError crcconfig([1, 0, 1, 0])
            @test_throws ArgumentError crcconfig([RRCFilters._CRC_MAX_POLYNOMIAL_DEGREE + 1, 0])

            @test_throws ArgumentError crcconfig([16, 12, 5, 0]; initial_conditions=2)
            @test_throws ArgumentError crcconfig([16, 12, 5, 0]; initial_conditions=[1, 0, 1])
            @test_throws ArgumentError crcconfig([16, 12, 5, 0]; initial_conditions=[1, 0, 1, 2])
            @test_throws ArgumentError crcconfig([16, 12, 5, 0]; final_xor=2)
            @test_throws ArgumentError crcconfig([16, 12, 5, 0]; final_xor=[1, 0])

            @test_throws ArgumentError crcconfig([16, 12, 5, 0]; checksums_per_frame=true)
            @test_throws ArgumentError crcconfig([16, 12, 5, 0]; checksums_per_frame=0)
            @test_throws ArgumentError crcconfig([16, 12, 5, 0]; checksums_per_frame=-1)
        end

        @testset "dual polynomial form parsing agreement" begin
            power_form = crcconfig([16, 12, 5, 0])
            binary_form = crcconfig([1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
            @test power_form == binary_form

            degree0 = crcconfig([1])
            @test degree0.num_bits == 0
            @test degree0.polynomial == Bool[true]
        end

        @testset "initial_conditions/final_xor accept Bool naturally (bit content, not a count)" begin
            @test crcconfig([16, 12, 5, 0]; initial_conditions=true) == crcconfig([16, 12, 5, 0]; initial_conditions=1)
            @test crcconfig([8, 7, 6, 4, 2, 0]; final_xor=BitVector([true, false, true, true, false, false, true, false])) ==
                  crcconfig([8, 7, 6, 4, 2, 0]; final_xor=[1, 0, 1, 1, 0, 0, 1, 0])
        end

        @testset "crcgenerate/crcdetect validation" begin
            cfg = crcconfig([8, 7, 6, 4, 2, 0])

            @test_throws ArgumentError crcgenerate([1, 0, 2, 1, 0, 0, 0, 0], cfg)
            @test_throws ArgumentError crcgenerate(bitrand(Xoshiro(1), 15), crcconfig([8, 7, 6, 4, 2, 0]; checksums_per_frame=2))
            @test_throws ArgumentError crcgenerate(bitrand(Xoshiro(2), 12), crcconfig([8, 7, 6, 4, 2, 0]; reflect_input_bytes=true))
            @test isempty(crcgenerate(Int[], cfg))

            @test_throws ArgumentError crcdetect(Int[], cfg)
            @test_throws ArgumentError crcdetect([1, 0, 2, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1], cfg)
            @test_throws ArgumentError crcdetect(bitrand(Xoshiro(3), 15), crcconfig([8, 7, 6, 4, 2, 0]; checksums_per_frame=2))
            @test_throws ArgumentError crcdetect(bitrand(Xoshiro(4), 8), cfg)
            @test_throws ArgumentError crcdetect(bitrand(Xoshiro(5), 20), crcconfig([8, 7, 6, 4, 2, 0]; reflect_input_bytes=true))
        end

        @testset "round trip across a config matrix" begin
            poly8 = [8, 7, 6, 4, 2, 0]
            poly16 = [16, 12, 5, 0]
            poly32 = [32, 26, 23, 22, 16, 12, 11, 10, 8, 7, 5, 4, 2, 1, 0]

            msg_a = bitrand(Xoshiro(201), 50)
            @test crc_round_trip_ok(msg_a, crcconfig(poly16))
            @test crc_round_trip_ok(msg_a, crcconfig(poly16; direct_method=true))
            @test crc_round_trip_ok(msg_a, crcconfig(poly16; initial_conditions=1))
            @test crc_round_trip_ok(msg_a, crcconfig(poly16; initial_conditions=1, direct_method=true))
            @test crc_round_trip_ok(msg_a, crcconfig(poly16; final_xor=1))
            @test crc_round_trip_ok(msg_a,
                crcconfig(poly16; initial_conditions=bitrand(Xoshiro(202), 16), final_xor=bitrand(Xoshiro(203), 16)))

            msg_b = bitrand(Xoshiro(204), 48)
            @test crc_round_trip_ok(msg_b, crcconfig(poly8; reflect_input_bytes=true, reflect_checksums=true))
            @test crc_round_trip_ok(msg_b,
                crcconfig(poly8; reflect_input_bytes=true, reflect_checksums=true, direct_method=true))

            msg_c = bitrand(Xoshiro(205), 40)
            @test crc_round_trip_ok(msg_c, crcconfig(poly8; checksums_per_frame=2))
            @test crc_round_trip_ok(msg_c, crcconfig(poly8; checksums_per_frame=4, direct_method=true))

            msg_d = bitrand(Xoshiro(206), 64)
            @test crc_round_trip_ok(msg_d,
                crcconfig(poly8; checksums_per_frame=2, reflect_input_bytes=true, reflect_checksums=true,
                    initial_conditions=1, final_xor=1))

            msg_e = bitrand(Xoshiro(207), 100)
            @test crc_round_trip_ok(msg_e, crcconfig(poly32))
            @test crc_round_trip_ok(msg_e, crcconfig(poly32; direct_method=true, initial_conditions=1, final_xor=1))

            @test crc_round_trip_ok(bitrand(Xoshiro(208), 20), crcconfig([1]))
        end

        @testset "direct vs indirect agreement at initial_conditions=0" begin
            msg = bitrand(Xoshiro(301), 64)
            for poly in ([8, 7, 6, 4, 2, 0], [16, 12, 5, 0], [16, 14, 1, 0])
                indirect_cfg = crcconfig(poly)
                direct_cfg = crcconfig(poly; direct_method=true)
                @test crcgenerate(msg, indirect_cfg) == crcgenerate(msg, direct_cfg)
            end

            diverge_indirect = crcconfig([16, 12, 5, 0]; initial_conditions=1)
            diverge_direct = crcconfig([16, 12, 5, 0]; initial_conditions=1, direct_method=true)
            @test crcgenerate(msg, diverge_indirect) != crcgenerate(msg, diverge_direct)
        end

        @testset "error detection" begin
            cfg = crcconfig([8, 7, 6, 4, 2, 0])
            msg = bitrand(Xoshiro(401), 32)
            codeword = crcgenerate(msg, cfg)
            for bit_index in eachindex(codeword)
                flipped = copy(codeword)
                flipped[bit_index] = !flipped[bit_index]
                _, err = crcdetect(flipped, cfg)
                @test any(err)
            end

            cfg2 = crcconfig([8, 7, 6, 4, 2, 0]; checksums_per_frame=2)
            msg2 = bitrand(Xoshiro(402), 32)
            codeword2 = crcgenerate(msg2, cfg2)
            frame_len = length(codeword2) ÷ 2

            flipped_first = copy(codeword2)
            flipped_first[3] = !flipped_first[3]
            _, err_first = crcdetect(flipped_first, cfg2)
            @test err_first == Bool[true, false]

            flipped_second = copy(codeword2)
            flipped_second[frame_len + 3] = !flipped_second[frame_len + 3]
            _, err_second = crcdetect(flipped_second, cfg2)
            @test err_second == Bool[false, true]
        end

        @testset "types, immutability, and generic indexing" begin
            cfg = crcconfig([16, 12, 5, 0])
            msg = bitrand(Xoshiro(501), 40)
            original_msg = copy(msg)

            codeword = @inferred crcgenerate(msg, cfg)
            @test msg == original_msg
            @test codeword isa BitVector

            original_codeword = copy(codeword)
            decoded_msg, err = @inferred crcdetect(codeword, cfg)
            @test codeword == original_codeword
            @test decoded_msg isa BitVector
            @test err isa BitVector
            @test decoded_msg == msg
            @test !any(err)

            zero_based_msg = ZeroBasedVector(collect(msg))
            zero_codeword = crcgenerate(zero_based_msg, cfg)
            @test zero_codeword == codeword

            zero_based_codeword = ZeroBasedVector(collect(codeword))
            zero_decoded_msg, zero_err = crcdetect(zero_based_codeword, cfg)
            @test zero_decoded_msg == msg
            @test !any(zero_err)
        end

        @test isempty(Test.detect_ambiguities(RRCFilters))
    end

    @testset "equalize" begin
        @testset "input validation" begin
            x = ComplexF64[1.0 + 0.5im, -0.3 + 0.8im, 0.0 + 0im, 0.9 - 0.2im, -0.6 - 0.6im]

            @test_throws ArgumentError equalize(x, fill(ComplexF64((1 + im) / sqrt(2)), length(x) + 5))

            @test_throws ArgumentError equalize(x; algorithm=:cma)
            @test_throws ArgumentError equalize(x; algorithm=:foo)
            @test (equalize(x; algorithm=:apa); true)

            @test_throws ArgumentError equalize(x; num_forward_taps=true)
            @test_throws ArgumentError equalize(x; num_forward_taps=0)
            @test_throws ArgumentError equalize(x; num_forward_taps=-1)

            @test_throws ArgumentError equalize(x; num_feedback_taps=true)
            @test_throws ArgumentError equalize(x; num_feedback_taps=-1)

            @test_throws ArgumentError equalize(x; reference_tap=true)
            @test_throws ArgumentError equalize(x; num_forward_taps=3, reference_tap=0)
            @test_throws ArgumentError equalize(x; num_forward_taps=3, reference_tap=4)

            @test_throws ArgumentError equalize(x; projection_order=true)
            @test_throws ArgumentError equalize(x; projection_order=0)
            @test_throws ArgumentError equalize(x; projection_order=-1)

            @test_throws ArgumentError equalize(x; step_size=true)
            @test_throws ArgumentError equalize(x; step_size=0.0)
            @test_throws ArgumentError equalize(x; step_size=-0.1)
            @test_throws ArgumentError equalize(x; step_size=NaN)
            @test_throws ArgumentError equalize(x; step_size=Inf)
            @test_throws ArgumentError equalize(x; algorithm=:apa, step_size=1.5)
            # step_size>1 stays legal for :lms/:rls -- the (0,1] bound is :apa-only.
            @test (equalize(x; algorithm=:lms, step_size=1.5); true)

            @test_throws ArgumentError equalize(x; forgetting_factor=true)
            @test_throws ArgumentError equalize(x; forgetting_factor=0.0)
            @test_throws ArgumentError equalize(x; forgetting_factor=1.0000001)
            @test_throws ArgumentError equalize(x; forgetting_factor=-0.5)
            @test_throws ArgumentError equalize(x; forgetting_factor=NaN)

            @test_throws ArgumentError equalize(x; initial_inverse_correlation=true)
            @test_throws ArgumentError equalize(x; initial_inverse_correlation=0.0)
            @test_throws ArgumentError equalize(x; initial_inverse_correlation=-1.0)
            @test_throws ArgumentError equalize(x; initial_inverse_correlation=NaN)

            @test_throws ArgumentError equalize(x; regularization=true)
            @test_throws ArgumentError equalize(x; regularization=0.0)
            @test_throws ArgumentError equalize(x; regularization=-1.0)
            @test_throws ArgumentError equalize(x; regularization=NaN)

            @test_throws ArgumentError equalize(x; constellation=ComplexF64[1.0 + 0im])
            @test_throws ArgumentError equalize(x; constellation=ComplexF64[])

            @test_throws ArgumentError equalize(ComplexF64[NaN + 0im, 1.0 + 0im])
            @test_throws ArgumentError equalize(ComplexF64[Inf + 0im, 1.0 + 0im])
        end

        @testset "num_feedback_taps=0 regression against pre-DFE behavior" begin
            rng = Xoshiro(91)
            x = randn(rng, ComplexF64, 256)
            training = randn(rng, ComplexF64, 64)
            for algorithm in (:lms, :rls, :apa)
                kwargs = algorithm === :apa ? (projection_order=3,) : NamedTuple()
                with_default = equalize(x, training; algorithm, num_forward_taps=7,
                                         reference_tap=3, kwargs...)
                with_explicit_zero = equalize(x, training; algorithm,
                                               num_forward_taps=7, num_feedback_taps=0,
                                               reference_tap=3, kwargs...)
                @test with_default == with_explicit_zero
            end
        end

        @testset "convergence / MSE-decreasing round trip" begin
            h = ComplexF64[2 * exp(im * 0.8 * pi), 3 * 0.25 * exp(-im * 0.4 * pi), 0.125 * exp(-im * 0.6 * pi)]
            for algorithm in (:lms, :rls)
                rng = Xoshiro(42)
                N = 2000
                bits = bitrand(rng, 2N)
                symbols = pskmod(bits, 4; phase_offset=pi / 4)
                channel_out = upfirdn(symbols, h, 1, 1)[1:N]
                noisy = awgn(rng, channel_out, 20.0)
                training = symbols[1:200]
                y, e, w = equalize(noisy, training; algorithm=algorithm, num_forward_taps=7, reference_tap=3)
                head_mse = mean(abs2, e[1:200])
                tail_mse = mean(abs2, e[(end - 200):end])
                @test tail_mse < head_mse / 2
            end
        end

        @testset "convergence / MSE-decreasing round trip — DFE" begin
            h = ComplexF64[2 * exp(im * 0.8 * pi), 3 * 0.25 * exp(-im * 0.4 * pi), 0.125 * exp(-im * 0.6 * pi)]
            for algorithm in (:lms, :rls)
                rng = Xoshiro(43)
                N = 2000
                bits = bitrand(rng, 2N)
                symbols = pskmod(bits, 4; phase_offset=pi / 4)
                channel_out = upfirdn(symbols, h, 1, 1)[1:N]
                noisy = awgn(rng, channel_out, 20.0)
                training = symbols[1:200]
                y, e, w = equalize(noisy, training; algorithm=algorithm, num_forward_taps=7, num_feedback_taps=3,
                                    reference_tap=3)
                @test length(w) == 10
                head_mse = mean(abs2, e[1:200])
                tail_mse = mean(abs2, e[(end - 200):end])
                @test tail_mse < head_mse / 2
            end
        end

        @testset "APA (:apa) — N=1 exactly reduces to NLMS" begin
            # Projection order one reduces exactly to NLMS, checked against
            # a standalone reference implementation.
            nlms_reference = (x, training, nf, reference_tap, mu, delta) -> begin
                T = ComplexF64
                n_train = length(training)
                delay = reference_tap - 1
                w = zeros(T, nf)
                u = zeros(T, nf)
                consts = RRCFilters._EQUALIZE_DEFAULT_CONSTELLATION
                y_out = Vector{T}(undef, length(x))
                e_out = Vector{T}(undef, length(x))
                for (step, sample) in enumerate(x)
                    for k in nf:-1:2
                        u[k] = u[k - 1]
                    end
                    u[1] = sample
                    y = sum(conj(w[k]) * u[k] for k in 1:nf)
                    y_out[step] = y
                    ref = if delay < step <= delay + n_train
                        T(training[step - delay])
                    else
                        best = consts[1]
                        bd = abs2(y - best)
                        for c in consts
                            d = abs2(y - c)
                            if d < bd
                                bd = d
                                best = c
                            end
                        end
                        best
                    end
                    e = ref - y
                    e_out[step] = e
                    norm2 = sum(abs2(u[k]) for k in 1:nf)
                    scale = mu / (delta + norm2)
                    for k in 1:nf
                        w[k] += scale * u[k] * conj(e)
                    end
                end
                return y_out, e_out, w
            end

            rng = Xoshiro(7)
            x = ComplexF64[randn(rng) + im * randn(rng) for _ in 1:500] ./ sqrt(2)
            training = ComplexF64[randn(rng) + im * randn(rng) for _ in 1:100] ./ sqrt(2)
            mu, delta, nf = 0.3, 1.0, 5

            y_ref, e_ref, w_ref = nlms_reference(x, training, nf, 3, mu, delta)
            y_apa, e_apa, w_apa = equalize(x, training; algorithm=:apa, num_forward_taps=nf, reference_tap=3,
                                            projection_order=1, step_size=mu, regularization=delta)

            @test y_apa ≈ y_ref atol=1e-10
            @test e_apa ≈ e_ref atol=1e-10
            @test w_apa ≈ w_ref atol=1e-10
        end

        @testset "convergence / MSE-decreasing round trip — APA" begin
            h = ComplexF64[2 * exp(im * 0.8 * pi), 3 * 0.25 * exp(-im * 0.4 * pi), 0.125 * exp(-im * 0.6 * pi)]
            rng = Xoshiro(44)
            N = 2000
            bits = bitrand(rng, 2N)
            symbols = pskmod(bits, 4; phase_offset=pi / 4)
            channel_out = upfirdn(symbols, h, 1, 1)[1:N]
            noisy = awgn(rng, channel_out, 20.0)
            training = symbols[1:200]
            y, e, w = equalize(noisy, training; algorithm=:apa, num_forward_taps=7, reference_tap=3,
                                projection_order=4, step_size=0.3)
            head_mse = mean(abs2, e[1:200])
            tail_mse = mean(abs2, e[(end - 200):end])
            @test tail_mse < head_mse / 2
        end

        @testset "convergence / MSE-decreasing round trip — APA with DFE" begin
            h = ComplexF64[2 * exp(im * 0.8 * pi), 3 * 0.25 * exp(-im * 0.4 * pi), 0.125 * exp(-im * 0.6 * pi)]
            rng = Xoshiro(45)
            N = 2000
            bits = bitrand(rng, 2N)
            symbols = pskmod(bits, 4; phase_offset=pi / 4)
            channel_out = upfirdn(symbols, h, 1, 1)[1:N]
            noisy = awgn(rng, channel_out, 20.0)
            training = symbols[1:200]
            y, e, w = equalize(noisy, training; algorithm=:apa, num_forward_taps=7, num_feedback_taps=3,
                                reference_tap=3, projection_order=4, step_size=0.3)
            @test length(w) == 10
            @test all(isfinite, y)
            head_mse = mean(abs2, e[1:200])
            tail_mse = mean(abs2, e[(end - 200):end])
            @test tail_mse < head_mse / 2
        end

        @testset "degenerate and edge cases" begin
            y, e, w = equalize(ComplexF64[])
            @test isempty(y) && isempty(e)
            @test length(w) == 5 && all(iszero, w)

            yr, er, wr = equalize(Float64[]; constellation=Float64[-1.0, 1.0])
            @test eltype(yr) === Float64
            @test length(wr) == 5 && all(iszero, wr)

            x1 = ComplexF64[1.0 + 0im, -1.0 + 0im, 0.9 + 0im, -0.8 + 0im, 1.0 + 0im, -0.9 + 0im]
            bipolar = ComplexF64[-1.0 + 0im, 1.0 + 0im]

            y1, e1, w1 = equalize(x1, ComplexF64[1.0 + 0im]; num_forward_taps=1, reference_tap=1, constellation=bipolar)
            @test length(w1) == 1
            @test all(isfinite, y1)

            y2, e2, w2 = equalize(x1; num_forward_taps=3, reference_tap=1, constellation=bipolar)
            @test all(isfinite, y2)
            @test length(w2) == 3

            y3, e3, w3 = equalize(x1, ComplexF64[1.0 + 0im]; num_forward_taps=4, reference_tap=4, constellation=bipolar)
            @test all(isfinite, y3)

            y4, e4, w4 = equalize(x1, ComplexF64[1.0 + 0im]; num_forward_taps=1, num_feedback_taps=1,
                                   reference_tap=1, constellation=bipolar)
            @test length(w4) == 2
            @test all(isfinite, y4)

            y5, e5, w5 = equalize(x1, ComplexF64[1.0 + 0im]; algorithm=:apa, num_forward_taps=1,
                                   reference_tap=1, constellation=bipolar, projection_order=1)
            @test length(w5) == 1
            @test all(isfinite, y5)
        end

        @testset "types, immutability, and generic indexing" begin
            x64 = ComplexF64[1.0 + 0.5im, -0.3 + 0.8im, 0.0 + 0im, 0.9 - 0.2im, -0.6 - 0.6im, 0.4 + 1.1im, 0.2 - 0.3im]
            ts64 = ComplexF64[(1 + im) / sqrt(2), (-1 + im) / sqrt(2)]
            original = copy(x64)

            y, e, w = equalize(x64, ts64; num_forward_taps=3, reference_tap=2)
            @test eltype(y) === ComplexF64
            @test eltype(e) === ComplexF64
            @test eltype(w) === ComplexF64

            x32 = ComplexF32.(x64)
            ts32 = ComplexF32.(ts64)
            constellation32 = ComplexF32.(RRCFilters._EQUALIZE_DEFAULT_CONSTELLATION)
            y32, e32, w32 = equalize(x32, ts32; num_forward_taps=3, reference_tap=2, constellation=constellation32)
            @test eltype(y32) === ComplexF32
            @test eltype(e32) === ComplexF32
            @test eltype(w32) === ComplexF32

            xr = Float64[1.0, -1.0, 0.9, -0.8, 1.0, -0.9, 0.8]
            tsr = Float64[1.0, -1.0]
            yr, er, wr = equalize(xr, tsr; num_forward_taps=3, reference_tap=2, constellation=Float64[-1.0, 1.0])
            @test eltype(yr) === Float64
            @test eltype(er) === Float64
            @test eltype(wr) === Float64

            equalize(x64, ts64; num_forward_taps=3, reference_tap=2)
            @test x64 == original

            zero_based_x = ZeroBasedVector(collect(x64))
            zero_based_ts = ZeroBasedVector(collect(ts64))
            @test equalize(zero_based_x, ts64; num_forward_taps=3, reference_tap=2) == equalize(x64, ts64; num_forward_taps=3, reference_tap=2)
            @test equalize(x64, zero_based_ts; num_forward_taps=3, reference_tap=2) == equalize(x64, ts64; num_forward_taps=3, reference_tap=2)

            @test @inferred(equalize(x64, ts64; num_forward_taps=3, reference_tap=2)) isa Tuple{Vector{ComplexF64},Vector{ComplexF64},Vector{ComplexF64}}

            y5, e5, w5 = equalize(x64, ts64; num_forward_taps=3, num_feedback_taps=2, reference_tap=2)
            @test eltype(y5) === ComplexF64
            @test length(w5) == 5
            @test @inferred(equalize(x64, ts64; num_forward_taps=3, num_feedback_taps=2, reference_tap=2)) isa
                  Tuple{Vector{ComplexF64},Vector{ComplexF64},Vector{ComplexF64}}

            y6, e6, w6 = equalize(x64, ts64; algorithm=:apa, num_forward_taps=3, reference_tap=2)
            @test eltype(y6) === ComplexF64
            @test eltype(e6) === ComplexF64
            @test eltype(w6) === ComplexF64
            @test @inferred(equalize(x64, ts64; algorithm=:apa, num_forward_taps=3, reference_tap=2)) isa
                  Tuple{Vector{ComplexF64},Vector{ComplexF64},Vector{ComplexF64}}

            y32_apa, e32_apa, w32_apa = equalize(x32, ts32; algorithm=:apa, num_forward_taps=3, reference_tap=2,
                                                  constellation=constellation32)
            @test eltype(y32_apa) === ComplexF32
            @test eltype(e32_apa) === ComplexF32
            @test eltype(w32_apa) === ComplexF32

            yr_apa, er_apa, wr_apa = equalize(xr, tsr; algorithm=:apa, num_forward_taps=3, reference_tap=2,
                                               constellation=Float64[-1.0, 1.0])
            @test eltype(yr_apa) === Float64
            @test eltype(er_apa) === Float64
            @test eltype(wr_apa) === Float64

            @test equalize(zero_based_x, ts64; algorithm=:apa, num_forward_taps=3, reference_tap=2) ==
                  equalize(x64, ts64; algorithm=:apa, num_forward_taps=3, reference_tap=2)
            @test equalize(x64, zero_based_ts; algorithm=:apa, num_forward_taps=3, reference_tap=2) ==
                  equalize(x64, ts64; algorithm=:apa, num_forward_taps=3, reference_tap=2)

            y7, e7, w7 = equalize(x64, ts64; algorithm=:apa, num_forward_taps=3, num_feedback_taps=2, reference_tap=2)
            @test eltype(y7) === ComplexF64
            @test length(w7) == 5
            @test all(isfinite, y7)
        end

        @test isempty(Test.detect_ambiguities(RRCFilters))
    end

    @testset "carrier-sync waveform link" begin
        function offset_link(bits, M, sps, span, beta, freq_offset, phase_offset, snr_db, modulation, rng)
            symbols = modulation == :qam ? qammod(bits, M) : pskmod(bits, M)
            txfilter = rcosdesign(beta, span, sps; shape=:sqrt)
            waveform = upfirdn(symbols, txfilter, sps)
            noisy_waveform = awgn(rng, waveform, snr_db)
            rxfilter = rcosdesign(beta, span, sps; shape=:sqrt)
            matched = upfirdn(noisy_waveform, rxfilter, 1, sps)
            symbol_count = length(symbols)
            recovered = matched[(span + 1):(span + symbol_count)]

            # Offset applied post-matched-filter (symbol rate), not on the
            # oversampled pre-filter waveform: a real carrier offset would
            # physically precede the receive matched filter, but doing so
            # here means the RRC filter (designed assuming a frequency-clean
            # signal) is no longer perfectly matched, adding its own
            # frequency-dependent distortion on top of what carriersync is
            # meant to demonstrate. This test's job is to show carriersync
            # integrates correctly with this library's own pulse-shaping/
            # channel/metrics functions -- not to model full RF-frontend
            # physical realism, which belongs in a dedicated channel model.
            offset_symbols = ComplexF64[s * cis(freq_offset * (i - 1) + phase_offset) for (i, s) in enumerate(recovered)]
            return offset_symbols, symbols
        end

        function best_demod_errors(corrected, M, bits, modulation, acquisition_symbols=0)
            bits_per_symbol = trailing_zeros(M)
            start = acquisition_symbols * bits_per_symbol + 1
            return minimum(0:(M - 1)) do rotation
                rotated = corrected .* cis(2pi * rotation / M)
                demod_bits = modulation == :qam ? qamdemod(rotated, M) : pskdemod(rotated, M)
                errors, _ = biterr(bits[start:end], demod_bits[start:end])
                return errors
            end
        end

        for (modulation, M, sps, span, beta, freq_offset, phase_offset, snr_db, seed) in
            ((:qam, 16, 4, 6, 0.25, 0.002, 0.4, 35.0, 501),)
            bits_per_symbol = trailing_zeros(M)
            bits = bitrand(Xoshiro(seed), 4096 * bits_per_symbol)

            recovered_symbols, ideal_symbols = offset_link(bits, M, sps, span, beta, freq_offset, phase_offset, snr_db, modulation, Xoshiro(seed + 1))
            corrected, phase_estimate = carriersync(recovered_symbols; modulation, sps=1)

            errors = best_demod_errors(corrected, M, bits, modulation, 300)
            @test errors == 0

            tail = (length(phase_estimate) - 199):length(phase_estimate)
            slope = (phase_estimate[tail[end]] - phase_estimate[tail[1]]) / (tail[end] - tail[1])
            @test isapprox(slope, freq_offset; atol=0.01)
        end

        let (modulation, M, sps, span, beta, freq_offset, phase_offset, snr_db, seed) = (:psk8, 8, 4, 6, 0.25, 0.0015, -0.3, 32.0, 502)
            bits_per_symbol = trailing_zeros(M)
            bits = bitrand(Xoshiro(seed), 4096 * bits_per_symbol)

            recovered_symbols, ideal_symbols = offset_link(bits, M, sps, span, beta, freq_offset, phase_offset, snr_db, modulation, Xoshiro(seed + 1))
            uncorrected_errors = best_demod_errors(recovered_symbols, M, bits, modulation)
            corrected, phase_estimate = carriersync(recovered_symbols; modulation, sps=1)
            corrected_errors = best_demod_errors(corrected, M, bits, modulation)
            @test corrected_errors == 0
            @test corrected_errors < uncorrected_errors

            tail = (length(phase_estimate) - 199):length(phase_estimate)
            slope = (phase_estimate[tail[end]] - phase_estimate[tail[1]]) / (tail[end] - tail[1])
            @test isapprox(slope, freq_offset; atol=0.01)
        end
    end
end
