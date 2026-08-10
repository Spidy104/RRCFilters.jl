using BenchmarkTools
using Random
using RRCFilters

const CONFIGURATIONS = (
    ("short", 0.25, 6, 4),
    ("typical", 0.35, 10, 8),
    ("long", 0.25, 32, 16),
)

function report_benchmark(label, beta, span, sps, shape; samples=1_000)
    rcosdesign(beta, span, sps; shape) # Compile before timing.
    trial = @benchmark rcosdesign($beta, $span, $sps; shape=$shape) samples=samples evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label), $(shape), $(span * sps + 1) taps: $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_awgn_benchmark(label, signal, snr_db; signal_power=nothing)
    rng = Xoshiro(1234)
    if isnothing(signal_power)
        awgn(rng, signal, snr_db) # Compile before timing.
        trial = @benchmark awgn($rng, $signal, $snr_db) samples=1_000 evals=1
    else
        awgn(rng, signal, snr_db; signal_power) # Compile before timing.
        trial = @benchmark awgn($rng, $signal, $snr_db; signal_power=$signal_power) samples=1_000 evals=1
    end
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_large_awgn_benchmark(label, signal, snr_db; signal_power=nothing)
    rng = Xoshiro(1234)
    if isnothing(signal_power)
        awgn(rng, signal, snr_db) # Compile before timing.
        trial = @benchmark awgn($rng, $signal, $snr_db) samples=100 evals=1
    else
        awgn(rng, signal, snr_db; signal_power) # Compile before timing.
        trial = @benchmark awgn($rng, $signal, $snr_db; signal_power=$signal_power) samples=100 evals=1
    end
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_qammod_benchmark(label, bits, M; gray=true)
    qammod(bits, M; gray) # Compile before timing.
    trial = @benchmark qammod($bits, $M; gray=$gray) samples=100 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_qamdemod_benchmark(label, symbols, M; gray=true)
    qamdemod(symbols, M; gray) # Compile before timing.
    trial = @benchmark qamdemod($symbols, $M; gray=$gray) samples=100 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_qamsoftdemod_benchmark(label, symbols, M, noise_variance; gray=true)
    qamsoftdemod(symbols, M; gray, noise_variance)
    trial = @benchmark qamsoftdemod($symbols, $M; gray=$gray, noise_variance=$noise_variance) samples=20 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function qam_awgn_link(rng, bits, M, snr_db; gray=true)
    symbols = qammod(bits, M; gray)
    noisy = awgn(rng, symbols, snr_db)
    return qamdemod(noisy, M; gray)
end

function report_qam_awgn_link_benchmark(label, bits, M, snr_db; gray=true)
    rng = Xoshiro(1234)
    qam_awgn_link(rng, bits, M, snr_db; gray) # Compile before timing.
    trial = @benchmark qam_awgn_link($rng, $bits, $M, $snr_db; gray=$gray) samples=20 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_pskmod_benchmark(label, bits, M; gray=true)
    pskmod(bits, M; gray) # Compile before timing.
    trial = @benchmark pskmod($bits, $M; gray=$gray) samples=100 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_pskdemod_benchmark(label, symbols, M; gray=true)
    pskdemod(symbols, M; gray) # Compile before timing.
    trial = @benchmark pskdemod($symbols, $M; gray=$gray) samples=100 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_psksoftdemod_benchmark(label, symbols, M, noise_variance; gray=true)
    psksoftdemod(symbols, M; gray, noise_variance)
    trial = @benchmark psksoftdemod($symbols, $M; gray=$gray, noise_variance=$noise_variance) samples=20 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function psk_awgn_link(rng, bits, M, snr_db; gray=true)
    symbols = pskmod(bits, M; gray)
    noisy = awgn(rng, symbols, snr_db)
    return pskdemod(noisy, M; gray)
end

function report_psk_awgn_link_benchmark(label, bits, M, snr_db; gray=true)
    rng = Xoshiro(1234)
    psk_awgn_link(rng, bits, M, snr_db; gray) # Compile before timing.
    trial = @benchmark psk_awgn_link($rng, $bits, $M, $snr_db; gray=$gray) samples=20 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_biterr_benchmark(label, tx_bits, rx_bits; mode=:bit, M=nothing)
    biterr(tx_bits, rx_bits; mode, M) # Compile before timing.
    trial = @benchmark biterr($tx_bits, $rx_bits; mode=$mode, M=$M) samples=1_000 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_evm_benchmark(label, ideal_symbols, rx_symbols; mode=:percent)
    evm(ideal_symbols, rx_symbols; mode) # Compile before timing.
    trial = @benchmark evm($ideal_symbols, $rx_symbols; mode=$mode) samples=1_000 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function qam_awgn_metrics_link(rng, bits, M, snr_db; gray=true)
    symbols = qammod(bits, M; gray)
    noisy = awgn(rng, symbols, snr_db)
    received_bits = qamdemod(noisy, M; gray)
    return (
        biterr(bits, received_bits),
        biterr(bits, received_bits; mode=:symbol, M),
        evm(symbols, noisy),
        evm(symbols, noisy; mode=:dB),
    )
end

function report_qam_awgn_metrics_link_benchmark(label, bits, M, snr_db; gray=true)
    rng = Xoshiro(1234)
    qam_awgn_metrics_link(rng, bits, M, snr_db; gray) # Compile before timing.
    trial = @benchmark qam_awgn_metrics_link($rng, $bits, $M, $snr_db; gray=$gray) samples=20 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_upfirdn_benchmark(label, x, h, p, q; samples=1_000)
    upfirdn(x, h, p, q) # Compile before timing.
    trial = @benchmark upfirdn($x, $h, $p, $q) samples=samples evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function waveform_link(rng, bits, M, sps, span, beta, snr_db)
    symbols = qammod(bits, M)
    txfilter = rcosdesign(beta, span, sps; shape=:sqrt)
    waveform = upfirdn(symbols, txfilter, sps)
    noisy_waveform = awgn(rng, waveform, snr_db)
    rxfilter = rcosdesign(beta, span, sps; shape=:sqrt)
    matched = upfirdn(noisy_waveform, rxfilter, 1, sps)
    symbol_count = length(symbols)
    recovered_symbols = matched[(span + 1):(span + symbol_count)]
    recovered_bits = qamdemod(recovered_symbols, M)
    return (
        biterr(bits, recovered_bits),
        evm(symbols, recovered_symbols),
    )
end

function report_waveform_link_benchmark(label, bits, M, sps, span, beta, snr_db)
    rng = Xoshiro(1234)
    waveform_link(rng, bits, M, sps, span, beta, snr_db) # Compile before timing.
    trial = @benchmark waveform_link($rng, $bits, $M, $sps, $span, $beta, $snr_db) samples=20 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_freqoffset_benchmark(label, x, freq_offset, phase_offset, sample_rate)
    freqoffset(x; freq_offset, phase_offset, sample_rate) # Compile before timing.
    trial = @benchmark freqoffset($x; freq_offset=$freq_offset, phase_offset=$phase_offset, sample_rate=$sample_rate) samples=100 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_carriersync_benchmark(label, x, modulation, sps)
    carriersync(x; modulation, sps) # Compile before timing.
    trial = @benchmark carriersync($x; modulation=$modulation, sps=$sps) samples=100 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_symbolsync_benchmark(label, x, sps, timing_error_detector; samples=100)
    symbolsync(x; sps, timing_error_detector) # Compile before timing.
    trial = @benchmark symbolsync($x; sps=$sps, timing_error_detector=$timing_error_detector) samples=samples evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_timingoffset_benchmark(label, x, delay)
    timingoffset(x; delay) # Compile before timing.
    trial = @benchmark timingoffset($x; delay=$delay) samples=100 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_agc_benchmark(label, x)
    agc(x) # Compile before timing.
    trial = @benchmark agc($x) samples=100 evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_poly2trellis_benchmark(label, constraint_length, code_generator; samples=1_000)
    poly2trellis(constraint_length, code_generator) # Compile before timing.
    trial = @benchmark poly2trellis($constraint_length, $code_generator) samples=samples evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_convenc_benchmark(label, bits, trellis; samples=100)
    convenc(bits, trellis) # Compile before timing.
    trial = @benchmark convenc($bits, $trellis) samples=samples evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_vitdec_benchmark(label, code, trellis, traceback_depth; mode=:trunc, decision_type=:hard,
                                  num_soft_bits=1, samples=10)
    vitdec(code, trellis, traceback_depth; mode, decision_type, num_soft_bits) # Compile before timing.
    trial = @benchmark vitdec($code, $trellis, $traceback_depth; mode=$mode, decision_type=$decision_type,
                               num_soft_bits=$num_soft_bits) samples=samples evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function fec_link(rng, bits, trellis, snr_db)
    coded, _ = convenc(bits, trellis)
    symbols = pskmod(coded, 2)
    noisy = awgn(rng, symbols, snr_db)
    received = pskdemod(noisy, 2)
    return vitdec(received, trellis, 20; mode=:trunc)
end

function report_fec_link_benchmark(label, bits, trellis, snr_db; samples=10)
    rng = Xoshiro(1234)
    fec_link(rng, bits, trellis, snr_db) # Compile before timing.
    trial = @benchmark fec_link($rng, $bits, $trellis, $snr_db) samples=samples evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_equalize_benchmark(label, x, training_symbols; algorithm=:lms, num_forward_taps=7,
                                    num_feedback_taps=0, reference_tap=3, projection_order=4,
                                    step_size=(algorithm === :apa ? 0.3 : 0.01), samples=10)
    equalize(x, training_symbols; algorithm, num_forward_taps, num_feedback_taps, reference_tap,
             projection_order, step_size) # Compile before timing.
    trial = @benchmark equalize($x, $training_symbols; algorithm=$algorithm, num_forward_taps=$num_forward_taps,
                                 num_feedback_taps=$num_feedback_taps, reference_tap=$reference_tap,
                                 projection_order=$projection_order, step_size=$step_size) samples=samples evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_crcgenerate_benchmark(label, msg, cfg; samples=10)
    crcgenerate(msg, cfg) # Compile before timing.
    trial = @benchmark crcgenerate($msg, $cfg) samples=samples evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

function report_crcdetect_benchmark(label, codeword, cfg; samples=10)
    crcdetect(codeword, cfg) # Compile before timing.
    trial = @benchmark crcdetect($codeword, $cfg) samples=samples evals=1
    best = minimum(trial)
    noise_ratio = maximum(trial).time / best.time
    println("$(label): $(round(best.time / 1_000; digits=2)) us (min; noise ratio $(round(noise_ratio; digits=2))x), $(best.memory) bytes, $(best.allocs) allocations")
end

println("Julia $(VERSION), $(Sys.CPU_NAME)")
println("Steady-state timings; compilation is excluded.")
for (label, beta, span, sps) in CONFIGURATIONS, shape in (:sqrt, :normal)
    report_benchmark(label, beta, span, sps, shape)
end
for shape in (:sqrt, :normal)
    report_benchmark("huge (1,000,001 taps)", 0.25, 1_000_000, 1, shape; samples=100)
end

println("AWGN steady-state timings; compilation and RNG construction are excluded.")
report_awgn_benchmark("AWGN real, 1,024 samples, measured power", ones(1_024), 10.0)
report_awgn_benchmark("AWGN real, 1,024 samples, supplied power", ones(1_024), 10.0; signal_power=1.0)
report_awgn_benchmark("AWGN complex, 1,024 samples, measured power", fill((1 + im) / sqrt(2), 1_024), 10.0)
report_large_awgn_benchmark("AWGN real, 1,000,000 samples, measured power", ones(1_000_000), 10.0)
report_large_awgn_benchmark("AWGN real, 1,000,000 samples, supplied power", ones(1_000_000), 10.0; signal_power=1.0)
report_large_awgn_benchmark("AWGN complex, 1,000,000 samples, measured power", fill((1 + im) / sqrt(2), 1_000_000), 10.0)

qam16_bits = bitrand(Xoshiro(16), 4 * 1_000_000)
qam64_bits = bitrand(Xoshiro(64), 6 * 1_000_000)
qam256_bits = bitrand(Xoshiro(256), 8 * 1_000_000)
qam16_symbols = qammod(qam16_bits, 16)
qam64_symbols = qammod(qam64_bits, 64)
qam256_symbols = qammod(qam256_bits, 256)

println("QAM steady-state timings; compilation and RNG construction are excluded.")
report_qammod_benchmark("16-QAM modulation, 1,000,000 symbols", qam16_bits, 16)
report_qamdemod_benchmark("16-QAM hard demodulation, 1,000,000 symbols", qam16_symbols, 16)
report_qammod_benchmark("64-QAM modulation, 1,000,000 symbols", qam64_bits, 64)
report_qamdemod_benchmark("64-QAM hard demodulation, 1,000,000 symbols", qam64_symbols, 64)
report_qammod_benchmark("256-QAM modulation, 1,000,000 symbols", qam256_bits, 256)
report_qamdemod_benchmark("256-QAM hard demodulation, 1,000,000 symbols", qam256_symbols, 256)
report_qamsoftdemod_benchmark("16-QAM max-log soft demodulation, 100,000 symbols",
                              qam16_symbols[1:100_000], 16, 0.1)
report_qamsoftdemod_benchmark("64-QAM max-log soft demodulation, 100,000 symbols",
                              qam64_symbols[1:100_000], 64, 0.1)

qam32_bits = bitrand(Xoshiro(32), 5 * 1_000_000)
qam128_bits = bitrand(Xoshiro(128), 7 * 1_000_000)
qam512_bits = bitrand(Xoshiro(512), 9 * 1_000_000)
qam32_symbols = qammod(qam32_bits, 32)
qam128_symbols = qammod(qam128_bits, 128)
qam512_symbols = qammod(qam512_bits, 512)

println("Cross-QAM steady-state timings; compilation and RNG construction are excluded.")
report_qammod_benchmark("32-QAM modulation, 1,000,000 symbols", qam32_bits, 32)
report_qamdemod_benchmark("32-QAM hard demodulation, 1,000,000 symbols", qam32_symbols, 32)
report_qammod_benchmark("128-QAM modulation, 1,000,000 symbols", qam128_bits, 128)
report_qamdemod_benchmark("128-QAM hard demodulation, 1,000,000 symbols", qam128_symbols, 128)
report_qammod_benchmark("512-QAM modulation, 1,000,000 symbols", qam512_bits, 512)
report_qamdemod_benchmark("512-QAM hard demodulation, 1,000,000 symbols", qam512_symbols, 512)

psk4_bits = bitrand(Xoshiro(4), 2 * 1_000_000)
psk8_bits = bitrand(Xoshiro(8), 3 * 1_000_000)
psk16_bits = bitrand(Xoshiro(160), 4 * 1_000_000)
psk4_symbols = pskmod(psk4_bits, 4)
psk8_symbols = pskmod(psk8_bits, 8)
psk16_symbols = pskmod(psk16_bits, 16)

println("PSK steady-state timings; compilation and RNG construction are excluded.")
report_pskmod_benchmark("QPSK modulation, 1,000,000 symbols", psk4_bits, 4)
report_pskdemod_benchmark("QPSK hard demodulation, 1,000,000 symbols", psk4_symbols, 4)
report_pskmod_benchmark("8-PSK modulation, 1,000,000 symbols", psk8_bits, 8)
report_pskdemod_benchmark("8-PSK hard demodulation, 1,000,000 symbols", psk8_symbols, 8)
report_pskmod_benchmark("16-PSK modulation, 1,000,000 symbols", psk16_bits, 16)
report_pskdemod_benchmark("16-PSK hard demodulation, 1,000,000 symbols", psk16_symbols, 16)
report_psksoftdemod_benchmark("QPSK max-log soft demodulation, 100,000 symbols",
                              psk4_symbols[1:100_000], 4, 0.1)
report_psksoftdemod_benchmark("8-PSK max-log soft demodulation, 100,000 symbols",
                              psk8_symbols[1:100_000], 8, 0.1)

large_psk_link_bits = bitrand(Xoshiro(1600), 4 * 100_000)
report_psk_awgn_link_benchmark("16-PSK/AWGN/16-PSK link, 100,000 symbols, 35 dB", large_psk_link_bits, 16, 35.0)

large_link64_bits = bitrand(Xoshiro(640), 6 * 100_000)
large_link256_bits = bitrand(Xoshiro(2560), 8 * 100_000)
report_qam_awgn_link_benchmark("64-QAM/AWGN/64-QAM link, 100,000 symbols, 40 dB", large_link64_bits, 64, 40.0)
report_qam_awgn_link_benchmark("256-QAM/AWGN/256-QAM link, 100,000 symbols, 45 dB", large_link256_bits, 256, 45.0)

metric_tx_bits = bitrand(Xoshiro(700), 1_000_000)
metric_rx_bits = copy(metric_tx_bits)
for index in 1:997:length(metric_rx_bits)
    metric_rx_bits[index] = !metric_rx_bits[index]
end
metric_evm_bits = bitrand(Xoshiro(701), 6 * 1_000_000)
metric_evm_symbols = qammod(metric_evm_bits, 64)
metric_evm_received = awgn(Xoshiro(702), metric_evm_symbols, 24.0)

println("Metric steady-state timings; scalar metrics must not allocate input-sized temporaries.")
report_biterr_benchmark("biterr bit mode, 1,000,000 bits", metric_tx_bits, metric_rx_bits)
report_biterr_benchmark("biterr symbol mode, 16-QAM width, 1,000,000 bits", metric_tx_bits, metric_rx_bits; mode=:symbol, M=16)
report_biterr_benchmark("biterr symbol mode, 256-QAM width, 1,000,000 bits", metric_tx_bits, metric_rx_bits; mode=:symbol, M=256)
report_evm_benchmark("EVM percent, 1,000,000 symbols", metric_evm_symbols, metric_evm_received)
report_evm_benchmark("EVM dB, 1,000,000 symbols", metric_evm_symbols, metric_evm_received; mode=:dB)
report_qam_awgn_metrics_link_benchmark("256-QAM/AWGN/metrics link, 100,000 symbols, 45 dB", large_link256_bits, 256, 45.0)

rrc_txfilter = rcosdesign(0.25, 6, 4; shape=:sqrt)
rrc_rxfilter = rcosdesign(0.25, 6, 4; shape=:sqrt)
upfirdn_tx_symbols = qammod(bitrand(Xoshiro(800), 4 * 1_000_000), 16)
upfirdn_tx_waveform = upfirdn(upfirdn_tx_symbols, rrc_txfilter, 4)

println("upfirdn steady-state timings; compilation is excluded.")
report_upfirdn_benchmark("upfirdn TX pulse shaping, 1,000,000 16-QAM symbols, sps=4", upfirdn_tx_symbols, rrc_txfilter, 4, 1; samples=100)
report_upfirdn_benchmark("upfirdn RX matched filtering, 4,000,025 samples, sps=4", upfirdn_tx_waveform, rrc_rxfilter, 1, 4; samples=100)

long_rrc_filter = rcosdesign(0.25, 32, 16; shape=:sqrt)
long_tx_waveform = upfirdn(upfirdn_tx_symbols, long_rrc_filter, 16)
report_upfirdn_benchmark("upfirdn TX pulse shaping, 1,000,000 16-QAM symbols, sps=16, 513-tap filter", upfirdn_tx_symbols, long_rrc_filter, 16, 1; samples=10)
report_upfirdn_benchmark("upfirdn RX matched filtering, 513-tap filter, sps=16", long_tx_waveform, long_rrc_filter, 1, 16; samples=10)

link_bits = bitrand(Xoshiro(801), 6 * 10_000)
report_waveform_link_benchmark("64-QAM waveform-level pulse-shape/AWGN/matched-filter/metrics link, 10,000 symbols, sps=4, 30 dB", link_bits, 64, 4, 6, 0.25, 30.0)

carriersync_qam_symbols = qammod(bitrand(Xoshiro(900), 4 * 1_000_000), 16)
carriersync_bpsk_symbols = pskmod(bitrand(Xoshiro(901), 1_000_000), 2)
carriersync_psk8_symbols = pskmod(bitrand(Xoshiro(902), 3 * 1_000_000), 8)

println("carriersync steady-state timings; compilation is excluded.")
report_carriersync_benchmark("carriersync :qam, 1,000,000 symbols", carriersync_qam_symbols, :qam, 1)
report_carriersync_benchmark("carriersync :bpsk, 1,000,000 symbols", carriersync_bpsk_symbols, :bpsk, 1)
report_carriersync_benchmark("carriersync :psk8, 1,000,000 symbols", carriersync_psk8_symbols, :psk8, 1)

freqoffset_complex_samples = carriersync_qam_symbols
freqoffset_real_samples = ones(1_000_000)

println("freqoffset steady-state timings; compilation is excluded.")
report_freqoffset_benchmark("freqoffset complex, 1,000,000 samples", freqoffset_complex_samples, 1000.0, 0.5, 1.0e6)
report_freqoffset_benchmark("freqoffset real, 1,000,000 samples", freqoffset_real_samples, 137.5, 0.3, 48000.0)

symbolsync_sps4_symbols = qammod(bitrand(Xoshiro(1000), 4 * 1_000_000), 16)
symbolsync_sps4_waveform = upfirdn(symbolsync_sps4_symbols, rcosdesign(0.25, 6, 4; shape=:sqrt), 4)
symbolsync_sps8_waveform = upfirdn(symbolsync_sps4_symbols, rcosdesign(0.25, 6, 8; shape=:sqrt), 8)

println("symbolsync steady-state timings; compilation is excluded.")
report_symbolsync_benchmark("symbolsync :zero_crossing, sps=4, 4,000,006 samples", symbolsync_sps4_waveform, 4, :zero_crossing)
report_symbolsync_benchmark("symbolsync :gardner, sps=4, 4,000,006 samples", symbolsync_sps4_waveform, 4, :gardner)
report_symbolsync_benchmark("symbolsync :early_late, sps=4, 4,000,006 samples", symbolsync_sps4_waveform, 4, :early_late)
report_symbolsync_benchmark("symbolsync :mueller_muller, sps=4, 4,000,006 samples", symbolsync_sps4_waveform, 4, :mueller_muller)
report_symbolsync_benchmark("symbolsync :gardner, sps=8, 8,000,006 samples", symbolsync_sps8_waveform, 8, :gardner; samples=20)

timingoffset_real_samples = ones(1_000_000)
timingoffset_complex_samples = fill((1 + im) / sqrt(2), 1_000_000)

println("timingoffset steady-state timings; compilation is excluded.")
report_timingoffset_benchmark("timingoffset real, 1,000,000 samples", timingoffset_real_samples, 1.7)
report_timingoffset_benchmark("timingoffset complex, 1,000,000 samples", timingoffset_complex_samples, 1.7)

agc_real_samples = ones(1_000_000)
agc_complex_samples = fill((1 + im) / sqrt(2), 1_000_000)

println("agc steady-state timings; compilation is excluded.")
report_agc_benchmark("agc real, 1,000,000 samples", agc_real_samples)
report_agc_benchmark("agc complex, 1,000,000 samples", agc_complex_samples)

fec_k7_trellis = poly2trellis(7, [0o171, 0o133])
fec_k9_trellis = poly2trellis(9, [0o561, 0o753])
fec_k7_bits = bitrand(Xoshiro(1100), 1_000_000)
fec_k9_bits = bitrand(Xoshiro(1101), 100_000)
fec_k7_code, _ = convenc(fec_k7_bits, fec_k7_trellis)
fec_k9_code, _ = convenc(fec_k9_bits, fec_k9_trellis)

println("FEC steady-state timings; compilation is excluded. vitdec's O(num_states x symbol_count) " *
        "predecessor-history memory is a disclosed first-version characteristic -- " *
        "the K=9 case is deliberately smaller-N than the " *
        "K=7 case for exactly this reason.")
report_poly2trellis_benchmark("poly2trellis K=7 (rate 1/2, 64 states)", 7, [0o171, 0o133])
report_poly2trellis_benchmark("poly2trellis K=9 (rate 1/2, 256 states)", 9, [0o561, 0o753])
report_convenc_benchmark("convenc K=7, 1,000,000 bits", fec_k7_bits, fec_k7_trellis)
report_convenc_benchmark("convenc K=9, 100,000 bits", fec_k9_bits, fec_k9_trellis)
report_vitdec_benchmark("vitdec K=7 (64 states), 1,000,000 symbols, :trunc", fec_k7_code, fec_k7_trellis, 20; samples=10)
report_vitdec_benchmark("vitdec K=9 (256 states), 100,000 symbols, :trunc", fec_k9_code, fec_k9_trellis, 20; samples=10)

fec_k7_soft_code = Int.(fec_k7_code)
fec_k7_unquant_code = Float64[b ? -1.0 : 1.0 for b in fec_k7_code]
fec_k7_llrs = 4 .* fec_k7_unquant_code
report_vitdec_benchmark("vitdec K=7 (64 states), 1,000,000 symbols, :trunc, decision_type=:soft", fec_k7_soft_code,
                         fec_k7_trellis, 20; mode=:trunc, decision_type=:soft, num_soft_bits=3, samples=10)
report_vitdec_benchmark("vitdec K=7 (64 states), 1,000,000 symbols, :trunc, decision_type=:unquant",
                         fec_k7_unquant_code, fec_k7_trellis, 20; mode=:trunc, decision_type=:unquant, samples=10)
report_vitdec_benchmark("vitdec K=7 (64 states), 1,000,000 symbols, :trunc, decision_type=:llr",
                         fec_k7_llrs, fec_k7_trellis, 20; mode=:trunc, decision_type=:llr, samples=10)

fec_link_bits = bitrand(Xoshiro(1102), 10_000)
report_fec_link_benchmark("FEC K=7 encode/BPSK/AWGN/demod/vitdec link, 10,000 bits, 2 dB", fec_link_bits, fec_k7_trellis, 2.0; samples=10)

equalize_bits = bitrand(Xoshiro(1300), 2_000_000)
equalize_symbols = pskmod(equalize_bits, 4; phase_offset=pi / 4)
equalize_h = ComplexF64[2 * exp(im * 0.8 * pi), 3 * 0.25 * exp(-im * 0.4 * pi), 0.125 * exp(-im * 0.6 * pi)]
equalize_channel_out = upfirdn(equalize_symbols, equalize_h, 1, 1)[1:1_000_000]
equalize_noisy = awgn(Xoshiro(1301), equalize_channel_out, 20.0)
equalize_training = equalize_symbols[1:1_000]

println("equalize steady-state timings; compilation is excluded. RLS's O((num_forward_taps+num_feedback_taps)^2)- " *
        "per-step cost (inverse-correlation matrix update) is a disclosed, deliberate first-version " *
        "characteristic -- kept at a modest tap count here for exactly this reason.")
report_equalize_benchmark("equalize :lms, NumForwardTaps=7, 1,000,000 symbols", equalize_noisy, equalize_training; algorithm=:lms, num_forward_taps=7, reference_tap=3, samples=10)
report_equalize_benchmark("equalize :rls, NumForwardTaps=7, 1,000,000 symbols", equalize_noisy, equalize_training; algorithm=:rls, num_forward_taps=7, reference_tap=3, samples=10)
report_equalize_benchmark("equalize :rls, NumForwardTaps=15, 1,000,000 symbols", equalize_noisy, equalize_training; algorithm=:rls, num_forward_taps=15, reference_tap=3, samples=5)
report_equalize_benchmark("equalize DFE :lms, NumForwardTaps=5, NumFeedbackTaps=3, 1,000,000 symbols", equalize_noisy, equalize_training; algorithm=:lms, num_forward_taps=5, num_feedback_taps=3, reference_tap=3, samples=10)
report_equalize_benchmark("equalize DFE :rls, NumForwardTaps=5, NumFeedbackTaps=3, 1,000,000 symbols", equalize_noisy, equalize_training; algorithm=:rls, num_forward_taps=5, num_feedback_taps=3, reference_tap=3, samples=10)
println("equalize :apa timings use the same warmed Julia harness as LMS and RLS. " *
        "Its O(projection_order^2 + projection_order*num_taps) per-step cost only beats RLS's " *
        "O(num_taps^2) when projection_order is chosen well below num_taps -- shown deliberately at " *
        "both a favorable and an unfavorable ratio below, not cherry-picked.")
report_equalize_benchmark("equalize :apa, NumForwardTaps=7, ProjectionOrder=2, 1,000,000 symbols", equalize_noisy, equalize_training; algorithm=:apa, num_forward_taps=7, reference_tap=3, projection_order=2, samples=10)
report_equalize_benchmark("equalize :apa, NumForwardTaps=7, ProjectionOrder=4, 1,000,000 symbols (ProjectionOrder too close to NumForwardTaps)", equalize_noisy, equalize_training; algorithm=:apa, num_forward_taps=7, reference_tap=3, projection_order=4, samples=10)
report_equalize_benchmark("equalize :apa, NumForwardTaps=15, ProjectionOrder=4, 1,000,000 symbols", equalize_noisy, equalize_training; algorithm=:apa, num_forward_taps=15, reference_tap=3, projection_order=4, samples=10)
report_equalize_benchmark("equalize DFE :apa, NumForwardTaps=5, NumFeedbackTaps=3, ProjectionOrder=4, 1,000,000 symbols", equalize_noisy, equalize_training; algorithm=:apa, num_forward_taps=5, num_feedback_taps=3, reference_tap=3, projection_order=4, samples=10)

crc32_cfg = crcconfig([32, 26, 23, 22, 16, 12, 11, 10, 8, 7, 5, 4, 2, 1, 0])
crc32_direct_cfg = crcconfig([32, 26, 23, 22, 16, 12, 11, 10, 8, 7, 5, 4, 2, 1, 0]; direct_method=true)
crc_msg = bitrand(Xoshiro(1400), 1_000_000)
crc_codeword = crcgenerate(crc_msg, crc32_cfg)

println("CRC steady-state timings; compilation is excluded.")
report_crcgenerate_benchmark("crcgenerate CRC-32, 1,000,000 bits, indirect", crc_msg, crc32_cfg)
report_crcgenerate_benchmark("crcgenerate CRC-32, 1,000,000 bits, direct", crc_msg, crc32_direct_cfg)
report_crcdetect_benchmark("crcdetect CRC-32, 1,000,000 bits, indirect", crc_codeword, crc32_cfg)
report_crcdetect_benchmark("crcdetect CRC-32, 1,000,000 bits, direct", crc_codeword, crc32_direct_cfg)
