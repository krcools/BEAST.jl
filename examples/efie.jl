using CompScienceMeshes
using BEAST

Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
X = raviartthomas(Γ)

κ = 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

@hilbertspace j
@hilbertspace k
efie = @discretise t[k,j]==e[k]  j∈X k∈X
u = gmres(efie)

include("utils/postproc.jl")
include("utils/plotresults.jl")
