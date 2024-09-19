using CompScienceMeshes
using BEAST

Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
# Γ = meshsphere(radius=1.0, h=0.1)
@show length(Γ)
X = raviartthomas(Γ)

κ, η = 1.0, 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

@hilbertspace j
@hilbertspace k
efie = @discretise t[k,j]==e[k]  j∈X k∈X
u, ch = BEAST.gmres_ch(efie; restart=1500)

include("utils/postproc.jl")
include("utils/plotresults.jl")

