using Documenter
using DocumenterCitations
using BEAST

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:alpha)

makedocs(;
    modules=[BEAST],
    authors="Kristof Cools and contributors",
    sitename="BEAST.jl",
    checkdocs=:none,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://krcools.github.io/BEAST.jl",
        edit_link="master",
        assets=String[],
        collapselevel=1,
        sidebar_sitename=true,
        size_threshold=1_000_000,
    ),
    plugins=[bib],
    pages=[
        "Introduction" => "index.md",
        "Manual" => Any[
            "General Usage"=>"manual/usage.md",
            "Custom excitations"=>"manual/customexc.md",
            "Setting the Quadrature Strategy" => "manual/quadstrat.md",
            "Custom Quadrature Rules" => "manual/quadrule.md",
            "Custom Operators" => "manual/customop.md",
            "Application Examples"=>Any[
                "Time-Harmonic"=>Any[
                    "EFIE"=>"manual/examplesTH/efie.md",
                    "MFIE"=>"manual/examplesTH/mfie.md",
                ],
                "Time-Domain"=>Any[
                    "EFIE"=>"manual/examplesTD/tdefie.md",
                ],
            ],
            "System of Equations and Bilinear Forms" => "manual/bilinear.md",
        ],
        "Tutorials" => Any[
            "tutorials/efie.md",
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
                "Plane Wave"=>"excitations/planewave.md",
                "Dipole"=>"excitations/dipole.md",
                "Monopole"=>"excitations/monopole.md",
                "Linear Potential"=>"excitations/linearpotential.md",
            ],
            "Composed Operator"=>Any[
                "Composed Operator"=>"composedoperator/composedoperator.md"
            ],
            "Quasi Helmholtz Projectors"=>Any[
                "Low Frequency Stable PMCHWT"=>"projectors/pmchwt_lf.md",
                "Dense Grid Stable PMCHWT"=>"projectors/pmchwt_theta.md",
            ],
        ],
        "Basis Functions" => Any[
            "Overview"=>"bases/overview.md",
            "Spatial"=>Any[
                "Raviart Thomas"=>"bases/raviartthomas.md",
                "Buffa Christiansen"=>"bases/buffachristiansen.md",
                "Brezzi-Douglas-Marini"=>"bases/brezzidouglasmarini.md",
                "Graglia-Wilton-Peterson"=>"bases/gragliawiltonpeterson.md",
            ],
            "Temporal"=>Any[],
        ],
        "Geometry Representations" => Any[
            "Flat"=>"geometry/flat.md",
            "Curvilinear"=>Any[]
        ],
        "____________________________________" => Any[],
        "Internals" => Any[
            "Overview"=>"internals/overview.md",
            "The Matrix Assemble Routine"=>"internals/assemble.md",
            "Quadrature"=>"internals/quadstrat.md",
            "Parametric Domain"=>Any[],
            "Multithreading"=>Any[],
        ],
        "Contributing" => Any["contributing.md"],
        "References" => Any["references.md"],
        "API Reference" => "apiref.md",
    ],
)

deploydocs(;
    repo="github.com/krcools/BEAST.jl.git",
    target="build",
    # push_preview=true,
    forcepush=false,
    #devbranch = "feature/docs",
)
