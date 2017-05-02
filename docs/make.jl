using Documenter, BEAST

makedocs(clean=false)
deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/krcools/BEAST.jl.git",
)
