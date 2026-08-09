# Performance

The current measured record is maintained in the root `benchmarks.md` file. It
includes environment, methodology, minimum time, allocated bytes, and
allocation count.

## Benchmark correctly

- Warm every method before timing.
- Keep compilation, data construction, and RNG construction outside timed code.
- Interpolate variables into `BenchmarkTools` expressions.
- Report allocations with time.
- Compare repeated runs on the same machine; sub-microsecond rows can be noisy.
- Pair every optimization with an independent numerical or property test.

## Known cost centers

- Cross-QAM boundary handling is more expensive than square-axis decisions.
- 8-PSK carrier recovery performs more per-sample detector arithmetic than QAM
  or BPSK recovery.
- Full-block Viterbi predecessor history is the largest memory cost.
- RLS cost grows quadratically with tap count.
- APA is attractive only when projection order is well below total taps.

The benchmark harness is intentionally not part of ordinary package tests. Run
it for performance-sensitive changes:

```powershell
julia --project=dev --startup-file=no dev/benchmark.jl
```
