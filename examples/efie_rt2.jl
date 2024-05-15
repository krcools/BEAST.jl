using CompScienceMeshes
using BEAST

#Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
Γ = meshsphere(radius=1.0, h=0.2)
# Γ = CompScienceMeshes.meshmobius(h=0.035)
X = BEAST.raviartthomas2(Γ)

κ, η = 1.0, 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
# E = -η/(im*κ)*BEAST.CurlCurlGreen(κ, ẑ, point(2,0,0))
e = (n × E) × n

@hilbertspace j
@hilbertspace k
efie = @discretise t[k,j]==e[k]  j∈X k∈X
u, ch = BEAST.gmres_ch(efie; restart=1500)

include("utils/postproc.jl")
include("utils/plotresults.jl")
