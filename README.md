# RRCFilters.jl

RRCFilters is a focused Julia package for building and testing digital
communications links. It covers pulse shaping, modulation, channel
impairments, synchronization, coding, integrity checks, metrics, gain control,
and adaptive equalization without nonstandard runtime dependencies.

## What it provides

| Area | Functions |
|---|---|
| Pulse shaping and multirate FIR | `rcosdesign`, `upfirdn` |
| QAM and PSK | `qammod`, `qamdemod`, `pskmod`, `pskdemod` |
| Noise and impairments | `awgn`, `freqoffset`, `timingoffset` |
| Receiver recovery | `agc`, `carriersync`, `symbolsync`, `equalize` |
| Coding and integrity | `poly2trellis`, `convenc`, `vitdec`, `crcconfig`, `crcgenerate`, `crcdetect` |
| Link metrics | `biterr`, `evm` |

## Installation

RRCFilters currently targets Julia 1.12 and is not yet registered. From a local
checkout:

```julia
using Pkg
Pkg.develop(path="path/to/RRCFilters")
```

## Quick start

```julia
using RRCFilters
using Random

rng = Xoshiro(42)
bits = bitrand(rng, 40_000)
transmitted = qammod(bits, 16)
received = awgn(rng, transmitted, 20.0)
decoded = qamdemod(received, 16)

errors, ber = biterr(bits, decoded)
evm_percent = evm(transmitted, received)
```

The functions operate on ordinary Julia vectors, return newly allocated
results, validate non-finite and unsafe configurations, and support seeded RNGs
for reproducible simulations.

## Documentation

Start with the [documentation map](docs/README.md), then use:

- [Getting started](docs/src/getting-started.md)
- [Capability guide](docs/src/capabilities.md)
- [Practical link examples](docs/src/examples.md)
- [Algorithms and limitations](docs/src/algorithms.md)
- [API reference](docs/src/api.md)
- [Performance record](benchmarks.md)

## Verification

```powershell
julia --project=. --startup-file=no -e "using Pkg; Pkg.test()"
julia --project=dev --startup-file=no dev/smoke.jl
julia --project=dev --startup-file=no dev/benchmark.jl
```

The current native suite contains more than 2.29 million checks. Development
and documentation tooling are isolated under `dev/`, `test/`, and `docs/`, so
package users install only the runtime standard-library dependencies.

See [CHANGELOG.md](CHANGELOG.md) for release history and [LICENSE](LICENSE) for
license terms.
