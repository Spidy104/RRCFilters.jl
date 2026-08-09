# Changelog

All notable changes to RRCFilters are documented here.

## Unreleased

- Added a complete Documenter manual with capability guides, practical links,
  limitations, performance methodology, validation guidance, and grouped API
  reference.
- Isolated documentation dependencies in the `docs` workspace.

## 0.1.0 - 2026-08-09

- Added raised-cosine design and efficient `upfirdn` pulse shaping.
- Added square, rectangular, and cross QAM plus power-of-two PSK modulation.
- Added AWGN, frequency/timing impairments, BER/SER, and EVM.
- Added AGC, carrier recovery, and four timing-error detectors.
- Added convolutional coding, hard/soft/unquantized Viterbi decoding, and CRC.
- Added LMS, RLS, APA, and decision-feedback equalization.
- Hardened finite extreme values, public trellis validation, bounded
  configuration work, and scale-safe distance and metric calculations.
- Added native high-precision, property, boundary, allocation, waveform, smoke,
  and performance verification.
