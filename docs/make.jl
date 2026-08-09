using Documenter
using RRCFilters

makedocs(
    sitename="RRCFilters.jl",
    authors="Dibyojyoti Bhattacharjee",
    modules=[RRCFilters],
    checkdocs=:exports,
    doctest=true,
    remotes=nothing,
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link=nothing,
        repolink=nothing,
    ),
    pages=[
        "Home" => "index.md",
        "Getting started" => "getting-started.md",
        "Capabilities" => "capabilities.md",
        "Practical examples" => "examples.md",
        "Receiver guides" => [
            "Waveform links" => "guides/waveform-links.md",
            "Synchronization" => "guides/synchronization.md",
            "Coding and integrity" => "guides/coding-and-integrity.md",
            "Adaptive equalization" => "guides/equalization.md",
        ],
        "Algorithms and limitations" => "algorithms.md",
        "Performance" => "performance.md",
        "Verification" => "validation.md",
        "API reference" => "api.md",
    ],
)
