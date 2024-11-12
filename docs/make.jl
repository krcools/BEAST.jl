using Documenter
using DocumenterCitations
using BEAST

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :alpha)

makedocs(;
    modules = [BEAST],
    authors = "Kristof Cools and contributors",
    sitename = "BEAST.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://krcools.github.io/BEAST.jl",
        edit_link = "master",
        assets = String[],
        collapselevel = 1,
        sidebar_sitename = true,
    ),
    plugins = [bib],
    pages = [
        "Introduction" => "index.md",
        "Manual" => Any[
            "General Usage"=>"manual/usage.md",
            "Application Examples"=>Any[
                "Time-Hamronic EFIE"=>"manual/efie.md",
                "Time-Domain EFIE"=>"manual/tdefie.md",
            ],
        ],
        "Operators & Excitations" => Any[
            "Overview"=>"operators/overview.md",
            "Local Operators"=>Any["Identiy"=>"operators/identity.md",],
            "Boundary Integral Operators"=>Any[
                "Helmholtz"=>"operators/helmholtz.md",
                "Maxwell Single Layer"=>Any[
                    "Time Harmonic"=>"operators/maxwellsinglelayer.md",
                    "Time Domain"=>"operators/maxwellsinglelayer_td.md",
                ],
                "Maxwell Double Layer"=>Any[
                    "Time Harmonic"=>"operators/maxwelldoublelayer.md",
                    "Time Domain"=>"operators/maxwelldoublelayer_td.md",
                ],
            ],
            "Volume Integral Operators"=>Any[
                "Maxwell Single Layer"=>"operators/maxwellsinglelayerVIE.md",
                "Maxwell Double Layer"=>Any[],
            ],
            "Excitations"=>Any[
                "Plane Wave"=>"operators/planewave.md",
                "Dipole"=>Any[],
                "Monopole"=>Any[],
                "Linear Potential" => Any[],
            ],
        ],
        "Bases" => Any[
            "Spatial"=>Any[
                "Raviart Thomas"=>"bases/raviartthomas.md",
                "Buffa Christiansen"=>"bases/buffachristiansen.md",
                "Brezzi-Douglas-Marini"=>"bases/brezzidouglasmarini.md",
            ],
            "Temporal"=>Any[],
        ],
        "Geometry Representations" => Any["Flat"=>Any[], "Curvilinear"=>Any[]],
        "____________________________________" => Any[],
        "Internals" => Any[
            "Overview"=>"internals/overview.md",
            "Assembly"=>"internals/assemble.md",
            "Quadrature"=>"internals/quadstrat.md",
            "Parametric Domain"=>Any[],
            "Multithreading"=>Any[],
        ],
        "Contributing" => "contributing.md",
        "References" => "references.md",
        "API Reference" => "apiref.md",
    ],
)

deploydocs(;
    repo = "github.com/krcools/BEAST.jl.git",
    target = "build",
    push_preview = true,
    forcepush = true,
)
