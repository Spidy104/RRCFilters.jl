# RRCFilters Performance Record

Last verified: 2026-08-09

Environment: Julia 1.12.6, `znver4` CPU target, Windows. The complete harness
finished in 200.5 seconds.

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
| `rcosdesign`, 1,000,001 taps, square-root | 24.899 ms | 8,000,135 | 3 |
| AWGN, real, 1,000,000 samples, measured power | 3.082 ms | 8,000,087 | 4 |
| AWGN, complex, 1,000,000 samples, measured power | 15.513 ms | 16,000,071 | 3 |
| `upfirdn` TX, 1,000,000 symbols, 25 taps, `sps=4` | 10.723 ms | 64,001,127 | 10 |
| `upfirdn` RX, 4,000,025 samples, 25 taps, `sps=4` | 8.774 ms | 16,000,599 | 7 |
| `upfirdn` TX, 1,000,000 symbols, 513 taps, `sps=16` | 80.154 ms | 256,016,821 | 12 |
| `upfirdn` RX, 513 taps, `sps=16` | 90.979 ms | 16,005,406 | 8 |

The million-tap filter remains linear in output size with only three
allocations. `upfirdn` memory is dominated by its output, not polyphase
metadata; the single-sample extreme-`p` regression separately enforces bounded
configuration overhead.

## Modulation and metrics

All modulation rows contain 1,000,000 symbols. Metric rows contain 1,000,000
bits or symbols.

| Case | Minimum | Bytes | Allocations |
|---|---:|---:|---:|
| 16-QAM modulation | 14.298 ms | 16,000,135 | 4 |
| 16-QAM demodulation | 25.209 ms | 500,199 | 5 |
| 256-QAM modulation | 18.484 ms | 16,000,135 | 4 |
| 256-QAM demodulation | 27.272 ms | 1,000,167 | 5 |
| 32-QAM modulation | 6.729 ms | 16,000,151 | 4 |
| 32-QAM demodulation | 43.366 ms | 625,207 | 5 |
| 512-QAM modulation | 10.493 ms | 16,000,151 | 4 |
| 512-QAM demodulation | 47.283 ms | 1,125,239 | 5 |
| QPSK modulation | 9.881 ms | 16,000,199 | 5 |
| QPSK demodulation | 23.608 ms | 250,151 | 4 |
| 8-PSK modulation | 10.783 ms | 16,000,263 | 5 |
| 8-PSK demodulation | 34.059 ms | 375,143 | 4 |
| 16-PSK modulation | 12.009 ms | 16,000,391 | 5 |
| 16-PSK demodulation | 38.160 ms | 500,135 | 4 |
| `biterr`, bit mode | 825.5 us | 0 | 0 |
| `biterr`, 256-QAM symbol mode | 1.605 ms | 0 | 0 |
| EVM, percent | 1.261 ms | 0 | 0 |
| EVM, dB | 1.266 ms | 0 | 0 |

Scalar metrics allocate no input-sized temporaries. Cross-QAM demodulation is
the slowest modulation decision path because cut-corner and boundary samples
need stable nearest-point handling.

## Synchronization and impairments

| Case | Input size | Minimum | Bytes | Allocations |
|---|---:|---:|---:|---:|
| `carriersync`, QAM | 1,000,000 | 10.477 ms | 24,000,206 | 8 |
| `carriersync`, BPSK | 1,000,000 | 10.133 ms | 24,000,206 | 8 |
| `carriersync`, 8-PSK | 1,000,000 | 29.327 ms | 24,000,206 | 8 |
| `freqoffset`, complex | 1,000,000 | 3.258 ms | 16,016,142 | 6 |
| `freqoffset`, real | 1,000,000 | 3.126 ms | 16,061,582 | 6 |
| `symbolsync`, zero crossing, `sps=4` | 4,000,006 | 91.686 ms | 49,600,718 | 12 |
| `symbolsync`, Gardner, `sps=4` | 4,000,006 | 90.194 ms | 49,600,718 | 12 |
| `symbolsync`, early-late, `sps=4` | 4,000,006 | 75.575 ms | 49,600,718 | 12 |
| `symbolsync`, Mueller-Muller, `sps=4` | 4,000,006 | 88.344 ms | 49,600,718 | 12 |
| `timingoffset`, real | 1,000,000 | 2.296 ms | 8,000,071 | 3 |
| `timingoffset`, complex | 1,000,000 | 8.614 ms | 16,000,071 | 3 |
| `agc`, real | 1,000,000 | 14.857 ms | 16,000,142 | 6 |
| `agc`, complex | 1,000,000 | 29.195 ms | 24,000,142 | 6 |

The 8-PSK carrier detector spends more arithmetic per sample than the QAM and
BPSK paths. That cost is accepted because the seeded 32 dB acquisition link
recovers with 0 errors across 11,388 post-acquisition bits.

## FEC, CRC, and equalization

| Case | Minimum | Bytes | Allocations |
|---|---:|---:|---:|
| `convenc`, K=7, 1,000,000 bits | 8.110 ms | 250,151 | 4 |
| `vitdec`, K=7 hard, 1,000,000 symbols | 381.269 ms | 264,126,565 | 16 |
| `vitdec`, K=7 soft, 1,000,000 symbols | 520.888 ms | 264,126,581 | 16 |
| `vitdec`, K=7 unquantized, 1,000,000 symbols | 398.922 ms | 264,126,581 | 16 |
| `crcgenerate`, CRC-32 indirect, 1,000,000 bits | 2.393 ms | 1,129,564 | 14 |
| `crcgenerate`, CRC-32 direct, 1,000,000 bits | 2.438 ms | 1,129,516 | 14 |
| `crcdetect`, CRC-32 indirect, 1,000,000 bits | 2.361 ms | 1,129,644 | 17 |
| `crcdetect`, CRC-32 direct, 1,000,000 bits | 2.187 ms | 1,129,660 | 17 |
| LMS, 7 forward taps, 1,000,000 symbols | 40.209 ms | 32,000,494 | 10 |
| RLS, 7 forward taps, 1,000,000 symbols | 288.136 ms | 32,001,886 | 18 |
| RLS, 15 forward taps, 1,000,000 symbols | 945.903 ms | 32,005,397 | 19 |
| APA, 7 taps, projection order 2 | 159.019 ms | 32,001,134 | 18 |
| APA, 7 taps, projection order 4 | 368.911 ms | 32,001,630 | 18 |
| APA, 15 taps, projection order 4 | 534.019 ms | 32,002,398 | 18 |
| DFE APA, 5 forward + 3 feedback taps, order 4 | 415.604 ms | 32,001,710 | 18 |

Viterbi's predecessor history is `O(num_states * symbol_count)` and is the
largest remaining memory cost. RLS is quadratic in total tap count. APA is
advantageous only when projection order stays well below tap count; the
7-tap/order-4 row is retained specifically to prevent cherry-picking.

## Verification paired with this snapshot

- Full package suite: 2,292,853/2,292,853 passed; test execution 1 minute
  50.1 seconds, 243.1 seconds total command wall time including environment
  setup and precompilation.
- Seeded end-to-end smoke suite: passed in 10.4 seconds.
- Complete benchmark harness: passed in 200.5 seconds.

Refresh this file only from a complete harness run on current code. Keep older
numbers out of the active record; version control is the history mechanism.
