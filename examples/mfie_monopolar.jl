using CompScienceMeshes, BEAST

Γ = readmesh(joinpath(dirname(@__FILE__),"sphere2.in"))
X = raviartthomas(Γ)
Y = buffachristiansen(Γ)
Z = raviartthomas(Γ, BEAST.Continuity{:none})

ϵ, μ, ω = 1.0, 1.0, 1.0; κ = ω * √(ϵ*μ)
K, N = Maxwell3D.doublelayer(wavenumber=κ), NCross()
NK, Id = BEAST.DoubleLayerRotatedMW3D(im*κ), Identity()
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
H = -1/(im*μ*ω)*curl(E)
h = (n × H) × n
q = n × H

@hilbertspace j
@hilbertspace m
@hilbertspace g

mfie1 = @discretise (K+0.5N)[m,j] == h[m]  j∈X m∈Y
mfie2 = @discretise (NK-0.5Id)[g,j] == q[m]  j∈Z g∈Z

u1 = gmres(mfie1)
u2 = gmres(mfie2)

fcr1, _ = facecurrents(u1, X)
fcr2, _ = facecurrents(u2, Z)

using LinearAlgebra
using Plots

plot(title="Compare current density")
plot!(norm.(fcr1), label="MxMIFE")
scatter!(norm.(fcr2), label="DGMFIE", c=:red)

include("utils/postproc.jl")
include("utils/plotresults.jl")
