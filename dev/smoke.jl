using RRCFilters
using Random
using Statistics

for beta in (0.0, 0.25, 1.0), shape in (:sqrt, :normal)
    coefficients = rcosdesign(beta, 6, 4; shape)
    @assert length(coefficients) == 25
    @assert coefficients == reverse(coefficients)
    @assert isapprox(sum(abs2, coefficients), 1.0; rtol=1e-14)
    @assert all(isfinite, coefficients)
end

sample_count = 10_000
real_signal = ones(sample_count)
real_noisy = awgn(Xoshiro(1), real_signal, 10.0)
@assert isapprox(10 * log10(mean(abs2, real_signal) / mean(abs2, real_noisy - real_signal)), 10.0; atol=0.3)

complex_signal = fill((1 + im) / sqrt(2), sample_count)
complex_noisy = awgn(Xoshiro(2), complex_signal, 8.0)
complex_noise = complex_noisy - complex_signal
@assert isapprox(10 * log10(mean(abs2, complex_signal) / mean(abs2, complex_noise)), 8.0; atol=0.3)
@assert isapprox(mean(abs2, real.(complex_noise)), mean(abs2, imag.(complex_noise)); rtol=0.1)

qam_bits = BitVector([false, false, false, true, true, false, true, true])
qam_symbols = qammod(qam_bits, 16)
@assert isapprox(mean(abs2, qam_symbols), 1.0; rtol=8eps(Float64))
@assert qamdemod(awgn(Xoshiro(3), qam_symbols, 30.0), 16) == qam_bits

qam256_bits = BitVector([false, false, false, true, true, false, true, true])
qam256_symbols = qammod(qam256_bits, 256)
@assert qamdemod(awgn(Xoshiro(4), qam256_symbols, 45.0), 256) == qam256_bits

function balanced_qam_bits(M, symbol_count)
    parameters = RRCFilters._qam_parameters(M)
    bits = falses(symbol_count * parameters.bits_per_symbol)
    for symbol_index in 1:symbol_count
        label = (symbol_index - 1) % M
        RRCFilters._integer_to_bits!(bits, (symbol_index - 1) * parameters.bits_per_symbol + 1, parameters.bits_per_symbol, label)
    end
    return bits
end

function smoke_qam_link_metrics(M, snr_db, seed)
    bits = balanced_qam_bits(M, 4_096)
    symbols = qammod(bits, M)
    noisy_symbols = awgn(Xoshiro(seed), symbols, snr_db)
    received_bits = qamdemod(noisy_symbols, M)
    bit_errors, ber = biterr(bits, received_bits)
    symbol_errors, ser = biterr(bits, received_bits; mode=:symbol, M)
    evm_percent = evm(symbols, noisy_symbols)
    evm_db = evm(symbols, noisy_symbols; mode=:dB)
    expected_evm_percent = 100 * 10^(-snr_db / 20)

    @assert mean(abs2, symbols) ≈ 1.0 rtol=64eps(Float64)
    @assert isapprox(evm_percent, expected_evm_percent; rtol=0.1)
    @assert 0.0 <= ber <= 1.0
    @assert 0.0 <= ser <= 1.0
    println("$(M)-QAM at $(snr_db) dB: BER=$(ber) ($(bit_errors) errors), SER=$(ser) ($(symbol_errors) errors), EVM=$(round(evm_percent; digits=3))% ($(round(evm_db; digits=3)) dB)")
end

smoke_qam_link_metrics(16, 20.0, 301)
smoke_qam_link_metrics(64, 24.0, 302)
smoke_qam_link_metrics(256, 28.0, 303)
smoke_qam_link_metrics(32, 22.0, 307)
smoke_qam_link_metrics(128, 26.0, 308)

psk_bits = BitVector([false, true, true, false])
psk_symbols = pskmod(psk_bits, 4)
@assert isapprox(mean(abs2, psk_symbols), 1.0; rtol=8eps(Float64))
@assert pskdemod(awgn(Xoshiro(5), psk_symbols, 20.0), 4) == psk_bits

psk8_bits = BitVector([false, true, true, false, true, false, false, true, true])
psk8_symbols = pskmod(psk8_bits, 8)
@assert pskdemod(awgn(Xoshiro(6), psk8_symbols, 30.0), 8) == psk8_bits

function balanced_psk_bits(M, symbol_count)
    parameters = RRCFilters._psk_parameters(M)
    bits = falses(symbol_count * parameters.bits_per_symbol)
    for symbol_index in 1:symbol_count
        label = (symbol_index - 1) % M
        RRCFilters._integer_to_bits!(bits, (symbol_index - 1) * parameters.bits_per_symbol + 1, parameters.bits_per_symbol, label)
    end
    return bits
end

function smoke_psk_link_metrics(M, snr_db, seed)
    bits = balanced_psk_bits(M, 4_096)
    symbols = pskmod(bits, M)
    noisy_symbols = awgn(Xoshiro(seed), symbols, snr_db)
    received_bits = pskdemod(noisy_symbols, M)
    bit_errors, ber = biterr(bits, received_bits)
    evm_percent = evm(symbols, noisy_symbols)

    @assert mean(abs2, symbols) ≈ 1.0 rtol=8eps(Float64)
    @assert 0.0 <= ber <= 1.0
    println("$(M)-PSK at $(snr_db) dB: BER=$(ber) ($(bit_errors) errors), EVM=$(round(evm_percent; digits=3))%")
end

smoke_psk_link_metrics(4, 18.0, 304)
smoke_psk_link_metrics(8, 24.0, 305)
smoke_psk_link_metrics(16, 30.0, 306)

function smoke_waveform_link(M, sps, span, beta, snr_db, seed)
    bits_per_symbol = trailing_zeros(M)
    bits = bitrand(Xoshiro(seed), 2_048 * bits_per_symbol)
    symbols = qammod(bits, M)
    txfilter = rcosdesign(beta, span, sps; shape=:sqrt)
    waveform = upfirdn(symbols, txfilter, sps)
    noisy_waveform = awgn(Xoshiro(seed + 1), waveform, snr_db)
    rxfilter = rcosdesign(beta, span, sps; shape=:sqrt)
    matched = upfirdn(noisy_waveform, rxfilter, 1, sps)
    recovered_symbols = matched[(span + 1):(span + length(symbols))]
    recovered_bits = qamdemod(recovered_symbols, M)

    bit_errors, ber = biterr(bits, recovered_bits)
    evm_percent = evm(symbols, recovered_symbols)
    @assert length(recovered_symbols) == length(symbols)
    println("Waveform-level $(M)-QAM link (sps=$sps, span=$span, beta=$beta) at $(snr_db) dB: " *
            "BER=$(ber) ($(bit_errors) errors), EVM=$(round(evm_percent; digits=3))%")
end

smoke_waveform_link(16, 4, 6, 0.25, 30.0, 501)

function best_rotation_errors(corrected, M, bits, demod)
    return minimum(0:(M - 1)) do rotation
        rotated = corrected .* cis(2pi * rotation / M)
        errors, _ = biterr(bits, demod(rotated, M))
        return errors
    end
end

function smoke_carriersync_link(modulation, M, freq_offset, phase_offset, snr_db, seed)
    bits_per_symbol = trailing_zeros(M)
    bits = bitrand(Xoshiro(seed), 4_096 * bits_per_symbol)
    symbols = modulation == :qam ? qammod(bits, M) : pskmod(bits, M)
    noisy_symbols = awgn(Xoshiro(seed + 1), symbols, snr_db)
    offset_symbols = freqoffset(noisy_symbols; freq_offset=freq_offset / (2pi), phase_offset)
    corrected, phase_estimate = carriersync(offset_symbols; modulation, sps=1)

    demod = modulation == :qam ? qamdemod : pskdemod
    tail = (length(phase_estimate) - 199):length(phase_estimate)
    slope = (phase_estimate[tail[end]] - phase_estimate[tail[1]]) / (tail[end] - tail[1])
    # Skips the initial PLL acquisition transient (the loop starts at
    # phase=0 and needs some symbols to pull in to a nonzero true phase
    # offset) before counting errors, matching real receiver practice.
    corrected_errors = best_rotation_errors(corrected[301:end], M, bits[(300 * bits_per_symbol + 1):end], demod)
    uncorrected_errors = best_rotation_errors(offset_symbols[301:end], M, bits[(300 * bits_per_symbol + 1):end], demod)
    @assert corrected_errors == 0

    println("carriersync $(modulation) M=$(M) at $(snr_db) dB: freq_offset=$(freq_offset) rad/sample, " *
            "tracked slope=$(round(slope; digits=6)), errors $(uncorrected_errors) -> $(corrected_errors) " *
            "(post-acquisition, out of $(length(bits) - 300 * bits_per_symbol) bits)")
end

smoke_carriersync_link(:qam, 16, 0.002, 0.4, 35.0, 601)
smoke_carriersync_link(:bpsk, 2, 0.001, -0.5, 20.0, 602)
smoke_carriersync_link(:psk8, 8, 0.0015, -0.3, 32.0, 603)

function smoke_symbolsync_link(ted, sps, span, beta, true_delay, snr_db, seed)
    bits = bitrand(Xoshiro(seed), 2 * 4_096)
    symbols = pskmod(bits, 4; phase_offset=pi / 4)
    txfilter = rcosdesign(beta, span, sps; shape=:sqrt)
    shaped = upfirdn(symbols, txfilter, sps)
    delayed = timingoffset(shaped; delay=true_delay)
    noisy = awgn(Xoshiro(seed + 1), delayed, snr_db)
    rxfilter = rcosdesign(beta, span, sps; shape=:sqrt)
    matched = upfirdn(noisy, rxfilter, 1, 1)

    recovered_symbols, timing_error = symbolsync(matched; sps, timing_error_detector=ted)
    tail = timing_error[(end-199):end]
    centered = mod.(tail .- tail[1] .+ 0.5, 1.0) .- 0.5
    lock_spread = maximum(centered) - minimum(centered)

    recovered_tail = recovered_symbols[(end-999):end]
    demodulated = pskdemod(recovered_tail, 4; phase_offset=pi / 4)
    remodulated = pskmod(demodulated, 4; phase_offset=pi / 4)
    evm_percent = evm(remodulated, recovered_tail)

    println("symbolsync $(ted) sps=$(sps) at $(snr_db) dB: true_delay=$(true_delay) samples, " *
            "lock spread=$(round(lock_spread; digits=6)), tail EVM=$(round(evm_percent; digits=3))%, " *
            "recovered $(length(recovered_symbols)) symbols (nominal $(length(matched)/sps))")
end

smoke_symbolsync_link(:zero_crossing, 4, 8, 0.25, 1.3, 30.0, 701)
smoke_symbolsync_link(:gardner, 4, 8, 0.25, -0.7, 30.0, 702)
smoke_symbolsync_link(:early_late, 4, 8, 0.25, 2.1, 30.0, 703)
smoke_symbolsync_link(:mueller_muller, 4, 8, 0.25, 0.4, 30.0, 704)

function smoke_agc_link(M, half, snr_db, low_scale, high_scale, seed)
    bits_per_symbol = trailing_zeros(M)
    symbol_count = 2 * half
    bits = bitrand(Xoshiro(seed), symbol_count * bits_per_symbol)
    symbols = qammod(bits, M)

    noisy = awgn(Xoshiro(seed + 1), symbols, snr_db; signal_power=1.0)
    swung = copy(noisy)
    swung[1:half] .*= low_scale
    swung[(half+1):end] .*= high_scale

    direct_demod = qamdemod(swung, M)
    _, direct_ber = biterr(bits, direct_demod)

    corrected, gain = agc(swung; desired_amplitude=1.0)
    corrected_demod = qamdemod(corrected, M)
    _, corrected_ber = biterr(bits, corrected_demod)

    println("agc M=$(M) at $(snr_db) dB: amplitude swing $(low_scale)x/$(high_scale)x, " *
            "BER without agc=$(round(direct_ber; digits=4)), with agc=$(round(corrected_ber; digits=4)), " *
            "gain range=[$(round(minimum(gain); digits=4)), $(round(maximum(gain); digits=4))]")
end

smoke_agc_link(16, 6_000, 25.0, 0.05, 20.0, 901)

function smoke_fec_link(K, gens, message_length, snr_db, seed)
    trellis = poly2trellis(K, gens)
    mem = trellis.constraint_length - 1
    bits = bitrand(Xoshiro(seed), message_length)
    bits_term = vcat(bits, falses(mem))

    coded, _ = convenc(bits_term, trellis)
    coded_symbols = pskmod(coded, 2)
    noisy_coded = awgn(Xoshiro(seed + 1), coded_symbols, snr_db)
    received_coded_bits = pskdemod(noisy_coded, 2)
    decoded = vitdec(received_coded_bits, trellis, 20; mode=:term)
    decoded_message = decoded[1:message_length]
    _, coded_ber = biterr(bits, decoded_message)

    real_samples = real.(noisy_coded)
    decoded_unquant = vitdec(real_samples, trellis, 20; mode=:term, decision_type=:unquant)
    _, coded_unquant_ber = biterr(bits, decoded_unquant[1:message_length])

    uncoded_symbols = pskmod(bits, 2)
    noisy_uncoded = awgn(Xoshiro(seed + 1), uncoded_symbols, snr_db)
    uncoded_demod = pskdemod(noisy_uncoded, 2)
    _, uncoded_ber = biterr(bits, uncoded_demod)

    println("FEC K=$(K) gens=$(gens) rate=1/$(trellis.n) BPSK at $(snr_db) dB: " *
            "uncoded BER=$(round(uncoded_ber; digits=4)), coded hard-decision BER=$(round(coded_ber; digits=4)), " *
            "coded unquant (soft) BER=$(round(coded_unquant_ber; digits=4))")
end

smoke_fec_link(7, [0o171, 0o133], 20_000, 2.0, 1001)
smoke_fec_link(7, [0o171, 0o133], 20_000, 0.0, 1002)

function smoke_equalize_link(algorithm, num_forward_taps, reference_tap, message_length, training_length, snr_db,
                              seed; num_feedback_taps=0, step_size=(algorithm === :apa ? 0.3 : 0.01), projection_order=4)
    h = ComplexF64[2 * exp(im * 0.8 * pi), 3 * 0.25 * exp(-im * 0.4 * pi), 0.125 * exp(-im * 0.6 * pi)]
    bits = bitrand(Xoshiro(seed), 2 * message_length)
    symbols = pskmod(bits, 4; phase_offset=pi / 4)
    channel_out = upfirdn(symbols, h, 1, 1)[1:message_length]
    noisy = awgn(Xoshiro(seed + 1), channel_out, snr_db)

    direct_demod = pskdemod(noisy, 4; phase_offset=pi / 4)
    _, direct_ber = biterr(bits, direct_demod)

    # equalize's output at step n approximates the transmitted symbol from
    # reference_tap-1 steps earlier (the equalizer's output latency,
    # unaffected by num_feedback_taps), so the demodulated bit stream must
    # be shifted by that many symbols before comparison, the same way this
    # project already trims RRC group delay (matched[(span+1):(end-span)])
    # or carriersync's acquisition transient before computing BER.
    delay = reference_tap - 1
    training = symbols[1:training_length]
    equalized, _, _ = equalize(noisy, training; algorithm=algorithm, num_forward_taps=num_forward_taps,
                                num_feedback_taps=num_feedback_taps, reference_tap=reference_tap,
                                step_size=step_size, projection_order=projection_order)
    corrected_demod = pskdemod(equalized, 4; phase_offset=pi / 4)
    aligned_bits = bits[1:(2 * (message_length - delay))]
    aligned_demod = corrected_demod[(2 * delay + 1):end]
    _, corrected_ber = biterr(aligned_bits, aligned_demod)

    println("equalize $(algorithm) NumForwardTaps=$(num_forward_taps) NumFeedbackTaps=$(num_feedback_taps) " *
            "ReferenceTap=$(reference_tap) QPSK multipath at $(snr_db) dB: " *
            "uncorrected BER=$(round(direct_ber; digits=4)), equalized BER=$(round(corrected_ber; digits=4))")
end

smoke_equalize_link(:lms, 7, 3, 20_000, 400, 15.0, 1101)
smoke_equalize_link(:rls, 7, 3, 20_000, 400, 15.0, 1102)
smoke_equalize_link(:lms, 5, 3, 20_000, 400, 15.0, 1103; num_feedback_taps=3)
smoke_equalize_link(:rls, 5, 3, 20_000, 400, 15.0, 1104; num_feedback_taps=3)
smoke_equalize_link(:apa, 7, 3, 20_000, 400, 15.0, 1105)
smoke_equalize_link(:apa, 5, 3, 20_000, 400, 15.0, 1106; num_feedback_taps=3)

function smoke_crc_link(K, gens, message_length, snr_db, seed)
    crc_cfg = crcconfig([16, 12, 5, 0])
    trellis = poly2trellis(K, gens)
    mem = trellis.constraint_length - 1

    bits = bitrand(Xoshiro(seed), message_length)
    protected = crcgenerate(bits, crc_cfg)
    protected_term = vcat(protected, falses(mem))

    coded, _ = convenc(protected_term, trellis)
    coded_symbols = pskmod(coded, 2)
    noisy_coded = awgn(Xoshiro(seed + 1), coded_symbols, snr_db)
    received_coded_bits = pskdemod(noisy_coded, 2)
    decoded = vitdec(received_coded_bits, trellis, 20; mode=:term)
    decoded_protected = decoded[1:length(protected)]

    _, post_viterbi_ber = biterr(protected, decoded_protected)
    recovered_msg, err = crcdetect(decoded_protected, crc_cfg)

    println("CRC+FEC link K=$(K) gens=$(gens) CRC-16 at $(snr_db) dB: " *
            "post-Viterbi BER=$(round(post_viterbi_ber; digits=5)), CRC flagged residual error=$(any(err)), " *
            "message bit-exact=$(recovered_msg == bits)")
end

smoke_crc_link(7, [0o171, 0o133], 2_000, 4.0, 1201)
smoke_crc_link(7, [0o171, 0o133], 2_000, -1.0, 1202)

println("RRCFilters smoke checks passed")
