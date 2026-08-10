# Capabilities

## Pulse shaping and multirate processing

`rcosdesign` designs unit-energy raised-cosine or root-raised-cosine FIR taps.
`upfirdn` performs zero-stuffed interpolation, full FIR convolution, and
decimation in one operation. Together they implement transmit pulse shaping
and receive matched filtering.

## Modulation

`qammod`, `qamdemod`, and `qamsoftdemod` support Gray or binary mapping for
square, rectangular, and corner-cut cross constellations. The soft demodulator
returns max-log bit LLRs.

`pskmod`, `pskdemod`, and `psksoftdemod` support power-of-two PSK orders, Gray
or binary labels, and an arbitrary finite phase offset. Soft-demodulation
orders are capped at 4096 because exact max-log partition searches scale with
constellation size.

## Channel impairments

`awgn` adds real or complex Gaussian noise at a requested SNR. It can measure
signal power or accept a supplied linear-power value.

`freqoffset` applies carrier frequency and phase rotation. `timingoffset`
introduces a deterministic fractional delay for recovery experiments.

## Receiver recovery

`agc` provides bounded log-domain automatic gain control. `carriersync`
supports QAM, BPSK, and axis-aligned 8-PSK phase detectors. `symbolsync`
supports zero-crossing, Gardner, early-late, and Mueller-Muller timing error
detectors.

`equalize` supports LMS, RLS, and affine-projection adaptation in linear or
decision-feedback configurations.

## Coding and integrity

`poly2trellis` builds rate-`1/n`, single-input, non-recursive convolutional
trellises. `convenc` encodes and `vitdec` performs hard, quantized-soft,
unquantized-sample, or LLR full-block Viterbi decoding in truncated or
terminated mode.

`crcconfig`, `crcgenerate`, and `crcdetect` support CRC degrees through 64,
direct or indirect operation, initial/final masks, byte reflection, checksum
reflection, and multiple checksums per frame.

## Metrics

`biterr` computes bit error rate or grouped symbol error rate. `evm` computes
normalized RMS error-vector magnitude in percent or decibels. Both use
scale-safe accumulation and the scalar result paths do not allocate
input-sized temporaries.
