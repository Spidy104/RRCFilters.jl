# [Coding and integrity](@id coding-guide)

## Convolutional coding

```julia
trellis = poly2trellis(7, [0o171, 0o133])
tail = falses(trellis.constraint_length - 1)
encoded, state = convenc(vcat(message, tail), trellis)
decoded = vitdec(encoded, trellis, 30; mode=:term)
```

The current trellis API is rate-`1/n`, binary, non-recursive, and
single-input-stream. Puncturing, erasures, recursive systematic codes, and
continuous streaming mode are outside scope.

Viterbi decision modes:

- `:hard`: exact zero/one observations
- `:soft`: integers from zero through `2^num_soft_bits - 1`
- `:unquant`: finite real matched-filter samples, with positive meaning zero
  and negative meaning one
- `:llr`: finite real log-likelihood ratios, with positive meaning zero and
  negative meaning one

Unquantized samples generally retain more information than thresholded hard
bits. LLRs from `qamsoftdemod` and `psksoftdemod` retain reliability across
higher-order constellations and connect directly to `decision_type=:llr`.
Terminated mode assumes the caller drove the encoder to state zero.

```julia
symbols = pskmod(encoded, 4; phase_offset=pi / 4)
received = awgn(rng, symbols, snr_db; signal_power=1.0)
llrs = psksoftdemod(received, 4; phase_offset=pi / 4,
                    noise_variance=10.0^(-snr_db / 10))
decoded = vitdec(llrs, trellis, 30; mode=:term, decision_type=:llr)
```

## CRC placement

CRC detects residual frame corruption; convolutional coding corrects likely
channel errors. A typical transmit order is:

```text
payload -> CRC append -> zero termination -> convolutional encode
```

The reverse receiver order is:

```text
Viterbi decode -> remove termination tail -> CRC detect -> payload
```

Do not treat a valid CRC as proof of zero errors; it is a finite-width detector
with a nonzero undetected-error probability.
