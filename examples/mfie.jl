using CompScienceMeshes, BEAST

Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
X, Y = raviartthomas(Γ), buffachristiansen(Γ)

ϵ, μ, ω = 1.0, 1.0, 1.0; κ, η = ω * √(ϵ*μ), √(μ/ϵ)
K, N = Maxwell3D.doublelayer(wavenumber=κ), NCross()
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
# E = -η/(im*κ)*BEAST.CurlCurlGreen(κ, ẑ, point(2,0,0))
H = -1/(im*μ*ω)*curl(E)
h = (n × H) × n

@hilbertspace j
@hilbertspace m
mfie = @discretise (K+0.5N)[m,j] == h[m]  j∈X m∈Y
u = gmres(mfie)

include("utils/postproc.jl")
include("utils/plotresults.jl")
