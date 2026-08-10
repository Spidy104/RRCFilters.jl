# Development Workflow

Instantiate the current development environment once:

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Run each verification path directly in its workspace:

```powershell
julia --project=. --startup-file=no -e "using Pkg; Pkg.test()"
julia --project=dev --startup-file=no dev/smoke.jl
julia --project=dev --startup-file=no dev/softcoded_link.jl
julia --project=dev --startup-file=no dev/benchmark.jl
julia --project=docs --startup-file=no docs/make.jl
```

`test/runtests.jl` contains unit, property, high-precision-reference,
extreme-numeric, allocation, inference, generic-indexing, and waveform-link tests.
`dev/smoke.jl` is the short seeded end-to-end gate. `dev/benchmark.jl` is the
complete warmed performance harness; it reports minimum time, timing noise,
allocated bytes, and allocation count while excluding compilation and RNG
construction from each timed call. `docs/make.jl` executes documentation
examples, checks exported API coverage, and renders the manual.

`dev/softcoded_link.jl` is the receiver-quality gate. It prints reproducible
CSV-formatted BER points for a pulse-shaped, timing- and carrier-impaired,
soft-decoded QPSK frame and fails if acquisition or coding-gain invariants do
not hold.

Always finish a change with at least the smoke command. Use the full test suite
for source changes and the complete benchmark harness for performance-sensitive
changes. Package precompilation and explicit warm-up inside `dev/benchmark.jl`
already keep compilation out of reported kernel timings.

`engineering-history.md` records durable design decisions, while
`benchmarks.md` is dated performance evidence.

The root `Project.toml` declares only package runtime dependencies. Test and
development dependencies live in `test/Project.toml`, `dev/Project.toml`, and
`docs/Project.toml`; the root workspace resolves them together into the single root manifest for
the active checkout without exposing them to package users.
