# Adaptive equalization

`equalize` returns `(equalized, error, weights)` and supports linear or
decision-feedback operation.

## Algorithms

- LMS: lowest computational cost; step size controls the convergence/noise
  tradeoff.
- RLS: faster convergence on difficult channels, but quadratic cost in total
  tap count.
- APA: projects over recent input vectors. It is useful when projection order
  is much smaller than tap count; projection order one is NLMS.

Set `num_feedback_taps=0` for a linear equalizer or a positive value for DFE.
Training symbols are used with `reference_tap - 1` samples of output latency;
after training, adaptation becomes decision-directed. An empty training vector
starts decision-directed operation immediately.

## Alignment

If `reference_tap=3`, compare the equalizer output starting two symbols later
with transmitted symbols ending two symbols earlier. Failure to compensate this
latency can make a converged equalizer appear broken.

## Tuning discipline

1. Verify channel/input normalization.
2. Start with a short supervised training segment.
3. Plot or compare head and tail MSE.
4. Confirm BER after latency alignment.
5. Add feedback taps only for postcursor ISI.
6. Increase RLS tap count or APA projection order only when measurements justify
   their cost.

The current implementation is symbol-rate and does not provide blind CMA,
fractionally spaced equalization, or persistent streaming state.
