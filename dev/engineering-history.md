# Engineering History

This is the compact record of decisions that still matter to the current code.
Obsolete experiments and comparison logs are intentionally omitted.

## Numerical design

- Raised-cosine filters use closed-form singularity limits, high-precision
  reference tests, symmetry checks, and unit-energy normalization.
- AWGN uses scale-safe power estimation and supports supplied power to avoid an
  unnecessary measurement pass.
- EVM, cross-QAM distance, and unquantized Viterbi metrics scale operands before
  squaring so every finite input has defined behavior.
- Frequency-offset rationalization is bounded; overflowed normalized frequency
  cannot enter an unbounded continued-fraction loop.

## Modulation and synchronization

- Square QAM uses direct axis quantization. Cross QAM uses fast rectangular bins
  where the nearest point is unambiguous and a stable distance search at cut
  corners and boundaries.
- PSK uses direct sector decisions for common orders. Empty modulation does not
  allocate a full constellation, even at very large valid orders.
- The carrier loop is a second-order proportional-integrator design. The 8-PSK
  phase detector was replaced with a direct decision-directed detector after
  waveform testing exposed poor acquisition quality.
- Symbol timing recovery uses a Farrow interpolator and selectable zero-crossing,
  Gardner, early-late, or Mueller-Muller detectors.

## Coding and adaptation

- Trellises are validated completely before any bounds-elided Viterbi kernel.
  Constraint length and generator count are capped before exponential shifts or
  allocations can occur.
- Hard, quantized-soft, and unquantized Viterbi decisions share survivor logic;
  one-bit soft decisions are exactly equivalent to hard decisions.
- CRC processing evolved from bit-serial state to packed-word and byte-table
  kernels while retaining independent randomized differential tests.
- The equalizer uses one combined forward/feedback tap vector for LMS, RLS, and
  affine projection. APA solves a small regularized Gram system in place;
  projection order one is checked against an independent NLMS implementation.

## Performance discipline

- Optimize only after a warmed, allocation-reporting benchmark identifies a
  real bottleneck.
- Keep compilation, setup, and RNG construction outside timed kernels.
- Validate every optimized path against an independent reference and boundary
  tests before accepting speedups.
- Avoid metadata or lookup tables sized by configuration when output size or
  work is much smaller.
