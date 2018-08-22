using CompScienceMeshes, BEAST

Γ = readmesh(joinpath(dirname(@__FILE__),"sphere2.in"))
X, Y = raviartthomas(Γ), buffachristiansen(Γ)

ϵ, μ, ω = 1.0, 1.0, 1.0; κ = ω * √(ϵ*μ)
NK, Id = BEAST.DoubleLayerRotatedMW3D(im*κ), Identity()
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
H = -1/(im*μ*ω)*curl(E)
h = (n × H) × n

@hilbertspace j
@hilbertspace m
mfie = @discretise (K+0.5N)[m,j] == h[m]  j∈X m∈Y
u = gmres(mfie)

include("utils/postproc.jl")
include("utils/plotresults.jl")
