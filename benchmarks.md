# RRCFilters Performance Record

Last verified: 2026-08-14

Environment: Julia 1.12.6, `znver4` CPU target, Windows. The complete harness
finished in 262.3 seconds.

Run it with:

```powershell
julia --project=dev --startup-file=no dev/benchmark.jl
```

## Method

All values are minimum steady-state times from `BenchmarkTools`. Package
compilation, input construction, and RNG construction are outside timed calls.
Every row reports allocated bytes and allocation count; input-sized output
arrays are included. A minimum is a regression anchor, not a latency guarantee.
Rows with high max-to-min noise should be compared only across repeated local
runs on the same machine.

## Filters, noise, and rate conversion

| Case | Minimum | Bytes | Allocations |
|---|---:|---:|---:|
| `rcosdesign`, 25 taps, square-root | 0.4 us | 256 | 2 |
| `rcosdesign`, 81 taps, square-root | 1.3 us | 704 | 2 |
| `rcosdesign`, 513 taps, square-root | 8.5 us | 4,231 | 3 |
| `rcosdesign`, 1,000,001 taps, square-root | 34.771 ms | 8,000,135 | 3 |
| AWGN, real, 1,000,000 samples, measured power | 3.771 ms | 8,000,087 | 4 |
| AWGN, complex, 1,000,000 samples, measured power | 21.584 ms | 16,000,071 | 3 |
| `upfirdn` TX, 1,000,000 symbols, 25 taps, `sps=4` | 15.653 ms | 64,001,127 | 10 |
| `upfirdn` RX, 4,000,025 samples, 25 taps, `sps=4` | 9.618 ms | 16,000,599 | 7 |
| `upfirdn` TX, 1,000,000 symbols, 513 taps, `sps=16` | 104.018 ms | 256,016,821 | 12 |
| `upfirdn` RX, 513 taps, `sps=16` | 135.659 ms | 16,005,406 | 8 |

The million-tap filter remains linear in output size with only three
allocations. `upfirdn` memory is dominated by its output, not polyphase
metadata; the single-sample extreme-`p` regression separately enforces bounded
configuration overhead.

## Modulation and metrics

Hard modulation rows contain 1,000,000 symbols. Soft-demodulation rows contain
100,000 symbols. Metric rows contain 1,000,000 bits or symbols.

| Case | Minimum | Bytes | Allocations |
|---|---:|---:|---:|
| 16-QAM modulation | 21.465 ms | 16,000,135 | 4 |
| 16-QAM hard demodulation | 36.058 ms | 500,199 | 5 |
| 16-QAM max-log soft demodulation, 100,000 symbols | 9.616 ms | 3,200,295 | 8 |
| 64-QAM max-log soft demodulation, 100,000 symbols | 18.708 ms | 4,800,295 | 8 |
| 256-QAM modulation | 27.146 ms | 16,000,135 | 4 |
| 256-QAM hard demodulation | 36.215 ms | 1,000,167 | 5 |
| 32-QAM modulation | 12.035 ms | 16,000,151 | 4 |
| 32-QAM hard demodulation | 40.398 ms | 625,207 | 5 |
| 512-QAM modulation | 19.415 ms | 16,000,151 | 4 |
| 512-QAM hard demodulation | 62.444 ms | 1,125,239 | 5 |
| QPSK modulation | 11.751 ms | 16,000,199 | 5 |
| QPSK hard demodulation | 33.020 ms | 250,151 | 4 |
| QPSK max-log soft demodulation, 100,000 symbols | 4.761 ms | 1,600,583 | 14 |
| 8-PSK modulation | 15.622 ms | 16,000,263 | 5 |
| 8-PSK hard demodulation | 45.391 ms | 375,143 | 4 |
| 8-PSK max-log soft demodulation, 100,000 symbols | 10.042 ms | 2,400,711 | 14 |
| 16-PSK modulation | 15.568 ms | 16,000,391 | 5 |
| 16-PSK hard demodulation | 49.492 ms | 500,135 | 4 |
| `biterr`, bit mode | 825.2 us | 0 | 0 |
| `biterr`, 256-QAM symbol mode | 1.616 ms | 0 | 0 |
| EVM, percent | 1.326 ms | 0 | 0 |
| EVM, dB | 1.354 ms | 0 | 0 |

Scalar metrics allocate no input-sized temporaries. Square-QAM soft metrics use
separable I/Q partition searches; cross-QAM and PSK use full constellation
partitions, with a 4096-order public cap to bound configuration work. Cross-QAM
hard demodulation is the slowest hard decision path because cut-corner and
boundary samples need stable nearest-point handling.

## Synchronization and impairments

| Case | Input size | Minimum | Bytes | Allocations |
|---|---:|---:|---:|---:|
| `carriersync`, QAM | 1,000,000 | 13.409 ms | 24,000,206 | 8 |
| `carriersync`, BPSK | 1,000,000 | 12.291 ms | 24,000,206 | 8 |
| `carriersync`, 8-PSK | 1,000,000 | 37.684 ms | 24,000,206 | 8 |
| `freqoffset`, complex | 1,000,000 | 3.447 ms | 16,016,142 | 6 |
| `freqoffset`, real | 1,000,000 | 3.718 ms | 16,061,582 | 6 |
| `symbolsync`, zero crossing, `sps=4` | 4,000,006 | 145.065 ms | 49,600,718 | 12 |
| `symbolsync`, Gardner, `sps=4` | 4,000,006 | 143.635 ms | 49,600,718 | 12 |
| `symbolsync`, early-late, `sps=4` | 4,000,006 | 116.845 ms | 49,600,718 | 12 |
| `symbolsync`, Mueller-Muller, `sps=4` | 4,000,006 | 130.329 ms | 49,600,718 | 12 |
| `timingoffset`, real | 1,000,000 | 2.369 ms | 8,000,071 | 3 |
| `timingoffset`, complex | 1,000,000 | 11.220 ms | 16,000,071 | 3 |
| `agc`, real | 1,000,000 | 19.476 ms | 16,000,142 | 6 |
| `agc`, complex | 1,000,000 | 40.148 ms | 24,000,142 | 6 |

The 8-PSK carrier detector spends more arithmetic per sample than the QAM and
BPSK paths. That cost is accepted because the seeded 32 dB acquisition link
recovers with 0 errors across 11,388 post-acquisition bits.

## FEC, CRC, and equalization

| Case | Minimum | Bytes | Allocations |
|---|---:|---:|---:|
| `convenc`, K=7, 1,000,000 bits | 9.931 ms | 250,151 | 4 |
| `vitdec`, K=7 hard, 1,000,000 symbols | 552.609 ms | 264,126,565 | 16 |
| `vitdec`, K=7 quantized-soft, 1,000,000 symbols | 460.429 ms | 264,126,597 | 16 |
| `vitdec`, K=7 unquantized, 1,000,000 symbols | 524.283 ms | 264,126,565 | 16 |
| `vitdec`, K=7 LLR, 1,000,000 symbols | 517.106 ms | 264,126,565 | 16 |
| `crcgenerate`, CRC-32 indirect, 1,000,000 bits | 3.614 ms | 1,129,548 | 14 |
| `crcgenerate`, CRC-32 direct, 1,000,000 bits | 3.412 ms | 1,129,564 | 14 |
| `crcdetect`, CRC-32 indirect, 1,000,000 bits | 2.990 ms | 1,129,628 | 17 |
| `crcdetect`, CRC-32 direct, 1,000,000 bits | 3.529 ms | 1,129,628 | 17 |
| LMS, 7 forward taps, 1,000,000 symbols | 57.035 ms | 32,000,494 | 10 |
| RLS, 7 forward taps, 1,000,000 symbols | 408.398 ms | 32,001,886 | 18 |
| RLS, 15 forward taps, 1,000,000 symbols | 1,377.116 ms | 32,005,397 | 19 |
| APA, 7 taps, projection order 2 | 237.649 ms | 32,001,134 | 18 |
| APA, 7 taps, projection order 4 | 558.787 ms | 32,001,630 | 18 |
| APA, 15 taps, projection order 4 | 775.608 ms | 32,002,398 | 18 |
| DFE APA, 5 forward + 3 feedback taps, order 4 | 562.304 ms | 32,001,710 | 18 |

Viterbi's predecessor history is `O(num_states * symbol_count)` and is the
largest remaining memory cost. LLR mode adds one block-wide reliability-scale
pass but no extra input-sized allocation. RLS is quadratic in total tap count.
APA is advantageous only when projection order stays well below tap count; the
7-tap/order-4 row is retained specifically to prevent cherry-picking.

## Verification paired with this snapshot

- Full package suite: 2,293,100/2,293,100 passed; test execution 2 minutes
  29.6 seconds, 361.9 seconds total command wall time including environment
  setup and precompilation.
- Seeded end-to-end smoke suite: passed in 13.7 seconds.
- Synchronized hard-versus-soft BER sweep and 40-trial stress gate: passed in
  9.0 seconds.
- Documentation build and executable examples: passed in 30.1 seconds.
- Complete benchmark harness: passed in 262.3 seconds.

Refresh this file only from a complete harness run on current code. Keep older
numbers out of the active record; version control is the history mechanism.
