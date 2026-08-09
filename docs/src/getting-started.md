# Getting started

## Requirements

RRCFilters 0.1 targets Julia 1.12. The package is currently unregistered, so
develop it from a checkout:

```julia
using Pkg
Pkg.develop(path="path/to/RRCFilters")
```

Then load the package normally:

```julia
using RRCFilters
```

## First reproducible link

```@example getting_started
using RRCFilters
using Random

rng = Xoshiro(42)
bits = bitrand(rng, 20_000)
tx = qammod(bits, 16)
rx = awgn(rng, tx, 20.0)
decoded = qamdemod(rx, 16)

errors, ber = biterr(bits, decoded)
(errors=errors, ber=ber, evm_percent=evm(tx, rx))
```

The explicit `rng` makes the result repeatable. `qammod` maps MSB-first groups
of four bits to unit-average-power 16-QAM symbols; `qamdemod` returns a
`BitVector` in the same ordering.

## Conventions

- Angles are radians.
- Frequencies and sample rates are hertz.
- SNR is decibels; supplied signal power is linear power.
- Modulators consume MSB-first bits.
- QAM constellations have unit average power; PSK constellations have unit
  modulus.
- Public operations allocate outputs and do not modify inputs.
- Boolean configuration values are rejected where an integer or real value is
  required.

## Choosing the next page

- Need a function inventory? Read [Capabilities](@ref).
- Need pulse shaping? Read [Waveform links](@ref).
- Need clock/carrier recovery? Read [Synchronization](@ref).
- Need FEC or CRC? Read [Coding and integrity guide](@ref coding-guide).
- Need channel correction? Read [Adaptive equalization](@ref).
