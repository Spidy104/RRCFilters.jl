# RRCFilters Documentation

The root README is the package overview. This directory contains the complete
manual.

## Reading order

1. [Getting started](src/getting-started.md)
2. [Capabilities](src/capabilities.md)
3. [Practical examples](src/examples.md)
4. Receiver guides:
   - [Waveform links](src/guides/waveform-links.md)
   - [Synchronization](src/guides/synchronization.md)
   - [Coding and integrity](src/guides/coding-and-integrity.md)
   - [Adaptive equalization](src/guides/equalization.md)
5. [Algorithms and limitations](src/algorithms.md)
6. [Performance](src/performance.md)
7. [Verification](src/validation.md)
8. [API reference](src/api.md)

## Build locally

From the package root:

```powershell
julia --project=docs --startup-file=no -e "using Pkg; Pkg.instantiate()"
julia --project=docs --startup-file=no docs/make.jl
```

Open `docs/build/index.html` after the build. Documenter is confined to the
documentation workspace and is not installed for package users.
