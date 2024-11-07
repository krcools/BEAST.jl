using Documenter
using DocumenterCitations
using BEAST

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:alpha)

makedocs(;
    modules=[BEAST],
    authors="Kristof Cools and contributors",
    sitename="BEAST.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://krcools.github.io/BEAST.jl",
        edit_link="master",
        assets=String[],
        collapselevel=1,
        sidebar_sitename=true,
    ),
    plugins=[bib],
    pages=[
        "Introduction" => "index.md",
        "Manual" => Any["General Usage" => "manual/usage.md", "Application Examples" => "manual/examples.md"],
        "Operators & Excitations" => Any[
            "Overview" => "operators/overview.md",
            "Local Operators" => Any["Identiy" => "operators/identity.md",],
            "Integral Operators" => Any[
                "Helmholtz" => "operators/helmholtz.md",
                "Maxwell Single Layer" => "operators/maxwellsinglelayer.md",
                "Maxwell Double Layer" => "operators/maxwelldoublelayer.md",
            ],
            "Excitations" => Any[
                "Plane Wave" => "operators/planewave.md",
            ],
        ],
        "Bases" => Any[
            "Spatial" => Any["Raviart Thomas" => "bases/raviartthomas.md", "Buffa Christiansen" => "bases/buffachristiansen.md"],
            "Temporal" => Any[],
        ],
        "Geometry Representations" => Any["Flat" => Any[], "Curvilinear" => Any[]],
        "____________________________________" => Any[],
        "Internals" => Any[
            "Multithreading" => Any[],
        ],
        "Contributing" => "contributing.md",
        "References" => "references.md",
        "API Reference" => "apiref.md",
    ],
)

deploydocs(; repo="github.com/krcools/BEAST.jl.git", target="build", push_preview=true, forcepush=true)
