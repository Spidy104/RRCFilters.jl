# Algorithms and limitations

## Numerical behavior

All public signal inputs must be finite. Kernels that square magnitudes use
scaled norms or distances so large finite floating-point values do not overflow
into invalid decisions. Configuration values used for shifts, table sizes, or
bounds-elided indexing are validated before entering optimized loops.

Concrete floating-point precision is preserved where supported. Integer and
rational signal inputs generally promote to `Float64`. Modulators return
`ComplexF64`; bit-producing functions return `BitVector`; soft demodulators
return `Vector{Float64}`.

Soft QAM and PSK demodulation uses the max-log approximation. For each bit it
finds the nearest constellation member in the zero and one partitions and
returns `(minimum_one_distance - minimum_zero_distance) / noise_variance`.
Positive LLRs favor zero. Ordinary inputs use a common-energy-free distance
score; extreme finite inputs fall back to high-precision arithmetic and clamp
only values that cannot be represented by `Float64`.

Square QAM separates into independent in-phase and quadrature searches, reducing
each symbol from a full constellation scan to two axis scans. Cross-QAM and PSK
retain the general partition search. PSK drops the common constant-energy term
exactly, preventing constellation-rounding bias for samples near zero.

## Stateless vector model

Every call starts with fresh loop/adaptation state. This makes experiments
reproducible and simple, but means the package does not yet expose persistent
streaming objects, chunk continuity, multichannel matrices, GPU kernels, or
in-place public APIs.

## Receiver limitations

- Soft demodulation is max-log rather than full log-MAP and is limited to
  constellation orders through 4096.
- Carrier loops require acquisition time and retain constellation phase
  ambiguity.
- `symbolsync` returns a variable number of symbols.
- `timingoffset` includes an intrinsic two-sample Farrow group delay.
- Full-block Viterbi stores predecessor history proportional to state count
  times symbol count; trellis construction is capped at constraint length 20.
- FEC does not include puncturing, erasures, continuous mode, or recursive
  systematic encoders.
- Equalization is symbol-rate and lacks blind CMA.

## Scope versus a full modem

RRCFilters supplies physical/link-layer building blocks, not framing policy.
Preambles, pilot design, packet detection, phase-ambiguity resolution,
interleaving, channel estimation, synchronization state machines, and protocol
retries remain application responsibilities.
