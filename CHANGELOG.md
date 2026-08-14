# Changelog

All notable changes to RRCFilters are documented here.

## Unreleased

- Hardened cross-QAM nearest-point decisions and PSK LLR signs at extreme
  finite magnitudes near overflow and underflow.
- Replaced the example receiver's fragile wide-pull-in decision loop with
  preamble-aided coarse and QPSK feedforward fine carrier recovery, backed by
  signed-offset, timing-delay, and multi-seed low-SNR gates.
- Reduced the maximum convolutional-code constraint length to 20 so public
  trellis construction cannot silently request hundreds of megabytes.
- Added max-log `qamsoftdemod` and `psksoftdemod` bit LLRs with finite-extreme
  hardening and bounded constellation-order work.
- Added `decision_type=:llr` Viterbi decoding with block-global metric scaling.
- Added a synchronized, pulse-shaped, soft-decoded QPSK receiver and seeded
  hard-versus-soft BER sweep with acquisition and quality gates.
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
