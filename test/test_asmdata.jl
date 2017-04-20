using CompScienceMeshes
using BEAST

meshfile = Pkg.dir("BEAST","test","sphere.in")
Γ = readmesh(meshfile)
X = raviartthomas(Γ)

els, ad = BEAST.assemblydata(X)

els2 = elements(geometry(X))
ad2 = assemblydata(X)

@assert ad2.data == ad.data
@assert length(els) == numcells(Γ)

resize!(X.fns, 1)

els, ad = BEAST.assemblydata(X)
@assert length(els) == 2
@assert size(ad.data) == (1,3,2)
