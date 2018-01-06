using CompScienceMeshes, BEAST
o, x, y, z = euclidianbasis(3)

Γ = readmesh(joinpath(dirname(@__FILE__),"sphere2.in"))
X = raviartthomas(Γ)

κ = 1.0
t = Maxwell3D.singlelayer(gamma=im*κ)
E = Maxwell3D.planewave(direction=z, polarization=x, wavenumber=κ)
e = (n × E) × n


@hilbertspace j
@hilbertspace k
efie = @discretise t[k,j]==e[k]  j∈X k∈X
u = gmres(efie)

include("utils/postproc.jl")
include("utils/plotresults.jl")
