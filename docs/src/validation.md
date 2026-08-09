# Verification

## Full package suite

```powershell
julia --project=. --startup-file=no -e "using Pkg; Pkg.test()"
```

The suite covers validation, fixed references, high-precision numerical
references, invariants, randomized differential checks, decision boundaries,
extreme finite values, type inference, immutability, generic indexing,
allocation regressions, and complete waveform links.

## End-to-end smoke suite

```powershell
julia --project=dev --startup-file=no dev/smoke.jl
```

This exercises QAM/PSK links, pulse shaping, carrier and timing acquisition,
AGC, convolutional coding, hard and unquantized Viterbi, CRC, and all equalizer
algorithms with deterministic seeds and quality assertions.

## Documentation

```powershell
julia --project=docs --startup-file=no docs/make.jl
```

The documentation build executes `@example` blocks and requires every exported
function to appear in the API reference.

## Performance

```powershell
julia --project=dev --startup-file=no dev/benchmark.jl
```

Do not replace correctness gates with benchmarks. A faster receiver that loses
lock quality or changes boundary decisions is a regression.
