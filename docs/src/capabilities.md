# Capabilities

## Pulse shaping and multirate processing

`rcosdesign` designs unit-energy raised-cosine or root-raised-cosine FIR taps.
`upfirdn` performs zero-stuffed interpolation, full FIR convolution, and
decimation in one operation. Together they implement transmit pulse shaping
and receive matched filtering.

## Modulation

`qammod` and `qamdemod` support Gray or binary mapping for orders 4, 8, 16, 32,
64, 128, 256, and 512. This includes square, rectangular, and corner-cut cross
constellations.

`pskmod` and `pskdemod` support power-of-two PSK orders, Gray or binary labels,
and an arbitrary finite phase offset. Demodulation is hard-decision; soft QAM
or PSK likelihoods are not part of the current API.

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
trellises. `convenc` encodes and `vitdec` performs hard, quantized-soft, or
unquantized full-block Viterbi decoding in truncated or terminated mode.

`crcconfig`, `crcgenerate`, and `crcdetect` support CRC degrees through 64,
direct or indirect operation, initial/final masks, byte reflection, checksum
reflection, and multiple checksums per frame.

## Metrics

`biterr` computes bit error rate or grouped symbol error rate. `evm` computes
normalized RMS error-vector magnitude in percent or decibels. Both use
scale-safe accumulation and the scalar result paths do not allocate
input-sized temporaries.
