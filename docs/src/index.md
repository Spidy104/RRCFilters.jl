# RRCFilters.jl

RRCFilters is a vector-oriented digital communications toolkit for Julia. Its
pieces are designed to compose into complete, reproducible baseband links while
remaining useful independently.

## Build a complete link

```text
bits -> CRC -> convolutional encoder -> modulation -> pulse shaping
     -> impairments/noise -> AGC -> matched filter -> timing recovery
     -> carrier recovery -> equalization -> demodulation -> Viterbi -> CRC
     -> BER/SER/EVM
```

You can use only the stages your experiment needs. Every public function is
stateless: inputs fully determine outputs, apart from the explicit RNG used by
`awgn`.

## Design goals

- Clear Julia vector APIs
- Reproducible seeded simulations
- Scale-safe finite-number behavior
- Bounded work for extreme valid configurations
- Independent numerical and property testing
- No nonstandard runtime dependencies

Use [Getting started](@ref) for the first link, [Capabilities](@ref) to choose
components, and [Algorithms and limitations](@ref) before treating a simulation
as a receiver design.
