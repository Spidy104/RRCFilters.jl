# Waveform links

A symbol-level simulation skips pulse shape, sample rate, and matched-filter
effects. Use a waveform link when timing recovery or occupied bandwidth matters.

## Transmit side

1. Generate bits.
2. Map them with `qammod` or `pskmod`.
3. Design root-raised-cosine taps with `rcosdesign`.
4. Interpolate and filter with `upfirdn(symbols, taps, sps, 1)`.

## Receive side

1. Apply AGC before loops if the amplitude range is uncertain.
2. Match-filter and decimate with `upfirdn(waveform, taps, 1, sps)` when timing
   is already aligned.
3. Remove the combined transmit/receive group delay: for equal RRC filters,
   take `matched[(span + 1):(span + symbol_count)]`.
4. Demodulate and compute BER/EVM.

If timing is unknown, keep the receive matched-filter output oversampled and
send it to `symbolsync` instead of immediately decimating it.

The executable `dev/softcoded_link.jl` example goes further: it pulse-shapes a
terminated rate-1/2 QPSK frame, adds timing, carrier, phase, and AWGN
impairments, performs Gardner timing recovery, uses a known acquisition word
for coarse carrier frequency and QPSK symmetry for fine frequency recovery,
acquires a sync word, resolves residual complex gain, produces max-log LLRs,
and Viterbi decodes the payload. Its seeded sweep prints CSV-formatted hard and
soft BER data.

## Avoiding misleading results

- Do not compare demodulated bits before removing group delay.
- Do not call symbol-rate SNR a waveform sample SNR without accounting for
  oversampling and filtering.
- Discard synchronization acquisition transients before BER.
- Use a seeded RNG and report the exact bit count; a zero-error short run is not
  a BER estimate.
