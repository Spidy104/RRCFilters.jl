using RRCFilters
using Random

function _softcoded_frame_start(received, synchronization)
    synchronization_length = length(synchronization)
    best_start = 0
    best_correlation = 0.0 + 0.0im
    best_score = -Inf
    synchronization_energy = sum(abs2, synchronization)

    for start in 1:(length(received) - synchronization_length + 1)
        correlation = 0.0 + 0.0im
        received_energy = 0.0
        for offset in 0:(synchronization_length - 1)
            sample = received[start + offset]
            correlation += sample * conj(synchronization[offset + 1])
            received_energy += abs2(sample)
        end
        score = abs2(correlation) / (received_energy * synchronization_energy)
        if score > best_score
            best_start = start
            best_correlation = correlation
            best_score = score
        end
    end

    gain = best_correlation / synchronization_energy
    return best_start, gain, best_score
end

function softcoded_link_trial(ebno_db::Real; message_length::Integer=20_000, seed::Integer=4100,
                              sps::Integer=4, span::Integer=8, beta::Real=0.25,
                              timing_delay::Real=0.65, frequency_offset::Real=0.0002,
                              phase_offset::Real=0.45)
    trellis = poly2trellis(7, [0o171, 0o133])
    tail_length = trellis.constraint_length - 1
    message = bitrand(Xoshiro(seed), message_length)
    terminated = vcat(message, falses(tail_length))
    coded, final_state = convenc(terminated, trellis)
    final_state == 0 || error("terminated encoder did not return to state zero")

    modulation_order = 4
    modulation_phase = pi / 4
    coded_symbols = pskmod(coded, modulation_order; phase_offset=modulation_phase)
    acquisition = pskmod(bitrand(Xoshiro(seed + 1), 2 * 512), modulation_order;
                         phase_offset=modulation_phase)
    synchronization = pskmod(bitrand(Xoshiro(seed + 2), 2 * 256), modulation_order;
                             phase_offset=modulation_phase)
    suffix = pskmod(bitrand(Xoshiro(seed + 3), 2 * 64), modulation_order;
                    phase_offset=modulation_phase)
    frame = vcat(acquisition, synchronization, coded_symbols, suffix)

    filter = rcosdesign(beta, span, sps; shape=:sqrt)
    waveform = upfirdn(frame, filter, sps)
    delayed = timingoffset(waveform; delay=timing_delay)
    offset = freqoffset(delayed; freq_offset=frequency_offset, phase_offset)

    code_rate = 1 / trellis.n
    information_bits_per_symbol = code_rate * trailing_zeros(modulation_order)
    noise_variance = inv(information_bits_per_symbol * 10.0^(Float64(ebno_db) / 10))
    waveform_snr_db = 10 * log10((1 / sps) / noise_variance)
    noisy = awgn(Xoshiro(seed + 4), offset, waveform_snr_db; signal_power=1 / sps)

    matched = upfirdn(noisy, filter)
    timed, _ = symbolsync(matched; sps, timing_error_detector=:gardner)
    corrected, _ = carriersync(timed; modulation=:qam, sps=1)
    synchronization_start, gain, synchronization_score =
        _softcoded_frame_start(corrected, synchronization)

    data_start = synchronization_start + length(synchronization)
    data_end = data_start + length(coded_symbols) - 1
    data_end <= length(corrected) || error("receiver did not recover the complete coded frame")
    normalized = corrected[data_start:data_end] ./ gain
    normalized_noise_variance = noise_variance / abs2(gain)

    hard_coded = pskdemod(normalized, modulation_order; phase_offset=modulation_phase)
    coded_errors, coded_ber = biterr(coded, hard_coded)
    hard_decoded = vitdec(hard_coded, trellis, 30; mode=:term)
    hard_errors, hard_ber = biterr(message, hard_decoded[1:message_length])

    llrs = psksoftdemod(normalized, modulation_order; phase_offset=modulation_phase,
                       noise_variance=normalized_noise_variance)
    soft_decoded = vitdec(llrs, trellis, 30; mode=:term, decision_type=:llr)
    soft_errors, soft_ber = biterr(message, soft_decoded[1:message_length])

    return (
        ebno_db=Float64(ebno_db),
        message_length=Int(message_length),
        hard_errors,
        hard_ber,
        coded_errors,
        coded_ber,
        data_evm=evm(coded_symbols, normalized),
        soft_errors,
        soft_ber,
        synchronization_score,
        synchronization_start,
        recovered_symbols=length(corrected),
    )
end

function softcoded_ber_sweep(ebno_values=(4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
                             message_length::Integer=20_000, seed::Integer=4100)
    results = [softcoded_link_trial(ebno; message_length, seed) for ebno in ebno_values]
    println("EbN0_dB,bits,hard_errors,hard_BER,soft_errors,soft_BER,sync_score")
    for result in results
        println("$(result.ebno_db),$(result.message_length),$(result.hard_errors),$(result.hard_ber)," *
                "$(result.soft_errors),$(result.soft_ber),$(result.synchronization_score)")
    end
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    results = softcoded_ber_sweep()
    all(result -> result.synchronization_score > 0.55, results) ||
        error("synchronization quality gate failed")
    all(result -> result.synchronization_start == results[1].synchronization_start, results) ||
        error("frame acquisition was not stable across the BER sweep")
    sum(result.soft_errors for result in results) < sum(result.hard_errors for result in results) ||
        error("soft decoding did not improve the aggregate BER sweep")
    all(results[index].soft_ber >= results[index + 1].soft_ber for index in 1:(length(results) - 1)) ||
        error("soft-decoded BER was not monotonic across the seeded sweep")
end
