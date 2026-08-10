# RRCFilters Performance Record

Last verified: 2026-08-10

Environment: Julia 1.12.6, `znver4` CPU target, Windows. The complete harness
finished in 222.8 seconds.

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
| `rcosdesign`, 513 taps, square-root | 13.4 us | 4,231 | 3 |
| `rcosdesign`, 1,000,001 taps, square-root | 37.564 ms | 8,000,135 | 3 |
| AWGN, real, 1,000,000 samples, measured power | 3.746 ms | 8,000,087 | 4 |
| AWGN, complex, 1,000,000 samples, measured power | 15.594 ms | 16,000,071 | 3 |
| `upfirdn` TX, 1,000,000 symbols, 25 taps, `sps=4` | 14.600 ms | 64,001,127 | 10 |
| `upfirdn` RX, 4,000,025 samples, 25 taps, `sps=4` | 10.388 ms | 16,000,599 | 7 |
| `upfirdn` TX, 1,000,000 symbols, 513 taps, `sps=16` | 101.532 ms | 256,016,821 | 12 |
| `upfirdn` RX, 513 taps, `sps=16` | 87.660 ms | 16,005,406 | 8 |

The million-tap filter remains linear in output size with only three
allocations. `upfirdn` memory is dominated by its output, not polyphase
metadata; the single-sample extreme-`p` regression separately enforces bounded
configuration overhead.

## Modulation and metrics

Hard modulation rows contain 1,000,000 symbols. Soft-demodulation rows contain
100,000 symbols. Metric rows contain 1,000,000 bits or symbols.

| Case | Minimum | Bytes | Allocations |
|---|---:|---:|---:|
| 16-QAM modulation | 14.629 ms | 16,000,135 | 4 |
| 16-QAM hard demodulation | 25.082 ms | 500,199 | 5 |
| 16-QAM max-log soft demodulation, 100,000 symbols | 5.910 ms | 3,200,295 | 8 |
| 64-QAM max-log soft demodulation, 100,000 symbols | 12.057 ms | 4,800,295 | 8 |
| 256-QAM modulation | 17.821 ms | 16,000,135 | 4 |
| 256-QAM hard demodulation | 27.147 ms | 1,000,167 | 5 |
| 32-QAM modulation | 6.909 ms | 16,000,151 | 4 |
| 32-QAM hard demodulation | 43.131 ms | 625,207 | 5 |
| 512-QAM modulation | 15.772 ms | 16,000,151 | 4 |
| 512-QAM hard demodulation | 73.869 ms | 1,125,239 | 5 |
| QPSK modulation | 11.885 ms | 16,000,199 | 5 |
| QPSK hard demodulation | 36.402 ms | 250,151 | 4 |
| QPSK max-log soft demodulation, 100,000 symbols | 4.920 ms | 1,600,583 | 14 |
| 8-PSK modulation | 14.990 ms | 16,000,263 | 5 |
| 8-PSK hard demodulation | 50.402 ms | 375,143 | 4 |
| 8-PSK max-log soft demodulation, 100,000 symbols | 9.103 ms | 2,400,711 | 14 |
| 16-PSK modulation | 16.548 ms | 16,000,391 | 5 |
| 16-PSK hard demodulation | 57.695 ms | 500,135 | 4 |
| `biterr`, bit mode | 847.0 us | 0 | 0 |
| `biterr`, 256-QAM symbol mode | 1.657 ms | 0 | 0 |
| EVM, percent | 1.379 ms | 0 | 0 |
| EVM, dB | 1.331 ms | 0 | 0 |

Scalar metrics allocate no input-sized temporaries. Square-QAM soft metrics use
separable I/Q partition searches; cross-QAM and PSK use full constellation
partitions, with a 4096-order public cap to bound configuration work. Cross-QAM
hard demodulation is the slowest hard decision path because cut-corner and
boundary samples need stable nearest-point handling.

## Synchronization and impairments

| Case | Input size | Minimum | Bytes | Allocations |
|---|---:|---:|---:|---:|
| `carriersync`, QAM | 1,000,000 | 10.081 ms | 24,000,206 | 8 |
| `carriersync`, BPSK | 1,000,000 | 9.718 ms | 24,000,206 | 8 |
| `carriersync`, 8-PSK | 1,000,000 | 29.311 ms | 24,000,206 | 8 |
| `freqoffset`, complex | 1,000,000 | 3.623 ms | 16,016,142 | 6 |
| `freqoffset`, real | 1,000,000 | 3.436 ms | 16,061,582 | 6 |
| `symbolsync`, zero crossing, `sps=4` | 4,000,006 | 149.413 ms | 49,600,718 | 12 |
| `symbolsync`, Gardner, `sps=4` | 4,000,006 | 132.132 ms | 49,600,718 | 12 |
| `symbolsync`, early-late, `sps=4` | 4,000,006 | 110.022 ms | 49,600,718 | 12 |
| `symbolsync`, Mueller-Muller, `sps=4` | 4,000,006 | 115.052 ms | 49,600,718 | 12 |
| `timingoffset`, real | 1,000,000 | 2.145 ms | 8,000,071 | 3 |
| `timingoffset`, complex | 1,000,000 | 8.697 ms | 16,000,071 | 3 |
| `agc`, real | 1,000,000 | 15.190 ms | 16,000,142 | 6 |
| `agc`, complex | 1,000,000 | 29.705 ms | 24,000,142 | 6 |

The 8-PSK carrier detector spends more arithmetic per sample than the QAM and
BPSK paths. That cost is accepted because the seeded 32 dB acquisition link
recovers with 0 errors across 11,388 post-acquisition bits.

## FEC, CRC, and equalization

| Case | Minimum | Bytes | Allocations |
|---|---:|---:|---:|
| `convenc`, K=7, 1,000,000 bits | 8.434 ms | 250,151 | 4 |
| `vitdec`, K=7 hard, 1,000,000 symbols | 395.622 ms | 264,126,565 | 16 |
| `vitdec`, K=7 quantized-soft, 1,000,000 symbols | 395.862 ms | 264,126,613 | 16 |
| `vitdec`, K=7 unquantized, 1,000,000 symbols | 403.511 ms | 264,126,581 | 16 |
| `vitdec`, K=7 LLR, 1,000,000 symbols | 396.766 ms | 264,126,581 | 16 |
| `crcgenerate`, CRC-32 indirect, 1,000,000 bits | 2.859 ms | 1,129,548 | 14 |
| `crcgenerate`, CRC-32 direct, 1,000,000 bits | 2.717 ms | 1,129,532 | 14 |
| `crcdetect`, CRC-32 indirect, 1,000,000 bits | 2.814 ms | 1,129,628 | 17 |
| `crcdetect`, CRC-32 direct, 1,000,000 bits | 2.960 ms | 1,129,660 | 17 |
| LMS, 7 forward taps, 1,000,000 symbols | 41.294 ms | 32,000,494 | 10 |
| RLS, 7 forward taps, 1,000,000 symbols | 281.885 ms | 32,001,886 | 18 |
| RLS, 15 forward taps, 1,000,000 symbols | 934.950 ms | 32,005,397 | 19 |
| APA, 7 taps, projection order 2 | 161.347 ms | 32,001,134 | 18 |
| APA, 7 taps, projection order 4 | 372.616 ms | 32,001,630 | 18 |
| APA, 15 taps, projection order 4 | 541.529 ms | 32,002,398 | 18 |
| DFE APA, 5 forward + 3 feedback taps, order 4 | 390.267 ms | 32,001,710 | 18 |

Viterbi's predecessor history is `O(num_states * symbol_count)` and is the
largest remaining memory cost. LLR mode adds one block-wide reliability-scale
pass but no extra input-sized allocation. RLS is quadratic in total tap count.
APA is advantageous only when projection order stays well below tap count; the
7-tap/order-4 row is retained specifically to prevent cherry-picking.

## Verification paired with this snapshot

- Full package suite: 2,292,998/2,292,998 passed; test execution 2 minutes
  22.3 seconds, 317.6 seconds total command wall time including environment
  setup and precompilation.
- Seeded end-to-end smoke suite: passed in 8.6 seconds.
- Synchronized hard-versus-soft BER sweep: passed in 4.1 seconds.
- Documentation build and executable examples: passed in 13.9 seconds.
- Complete benchmark harness: passed in 222.8 seconds.

Refresh this file only from a complete harness run on current code. Keep older
numbers out of the active record; version control is the history mechanism.
