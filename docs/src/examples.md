# Practical examples

## Pulse-shaped 16-QAM link

```@example waveform_example
using RRCFilters
using Random

rng = Xoshiro(501)
M, sps, span = 16, 4, 6
bits = bitrand(rng, 2_048 * trailing_zeros(M))
symbols = qammod(bits, M)
rrc = rcosdesign(0.25, span, sps)

waveform = upfirdn(symbols, rrc, sps, 1)
noisy = awgn(rng, waveform, 30.0)
matched = upfirdn(noisy, rrc, 1, sps)
recovered = matched[(span + 1):(span + length(symbols))]
decoded = qamdemod(recovered, M)

errors, ber = biterr(bits, decoded)
(errors=errors, ber=ber, evm_percent=evm(symbols, recovered))
```

The two root-raised-cosine filters produce `span` symbols of combined delay,
which is why symbol extraction begins at `span + 1`.

## Terminated convolutional code

```@example fec_example
using RRCFilters
using Random

trellis = poly2trellis(7, [0o171, 0o133])
message = bitrand(Xoshiro(77), 1_000)
terminated = vcat(message, falses(trellis.constraint_length - 1))
code, final_state = convenc(terminated, trellis)
decoded = vitdec(code, trellis, 30; mode=:term)

(final_state=final_state, exact=decoded == terminated)
```

Termination is explicit: append the zero tail before encoding and remove it
from the decoded message if the application does not want it.

## Soft-decoded QPSK

```@example soft_fec_example
using RRCFilters
using Random

trellis = poly2trellis(7, [0o171, 0o133])
message = bitrand(Xoshiro(91), 5_000)
terminated = vcat(message, falses(trellis.constraint_length - 1))
code, _ = convenc(terminated, trellis)
symbols = pskmod(code, 4; phase_offset=pi / 4)

snr_db = 2.0
received = awgn(Xoshiro(92), symbols, snr_db; signal_power=1.0)
hard_code = pskdemod(received, 4; phase_offset=pi / 4)
hard = vitdec(hard_code, trellis, 30; mode=:term)

llrs = psksoftdemod(received, 4; phase_offset=pi / 4,
                    noise_variance=10.0^(-snr_db / 10))
soft = vitdec(llrs, trellis, 30; mode=:term, decision_type=:llr)

hard_errors, hard_ber = biterr(message, hard[1:length(message)])
soft_errors, soft_ber = biterr(message, soft[1:length(message)])
(hard_errors=hard_errors, hard_ber=hard_ber,
 soft_errors=soft_errors, soft_ber=soft_ber)
```

`noise_variance` is the complex noise power at the demodulator input. The
soft result preserves bit reliability instead of discarding it at a decision
boundary.

## CRC-protected payload

```@example crc_example
using RRCFilters
using Random

cfg = crcconfig([16, 12, 5, 0])
payload = bitrand(Xoshiro(88), 512)
protected = crcgenerate(payload, cfg)
recovered, error = crcdetect(protected, cfg)

(exact=recovered == payload, checksum_error=any(error))
```

CRC detects residual corruption; it does not correct it. Place CRC around the
payload and FEC around the protected frame when combining both stages.
