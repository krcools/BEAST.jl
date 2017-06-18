using CompScienceMeshes, BEAST
o, x, y, z = euclidianbasis(3)

Γ = readmesh(joinpath(dirname(@__FILE__),"sphere2.in"))
X, Y = raviartthomas(Γ), buffachristiansen(Γ)

ω = κ = 1.0
μ = ϵ = 1.0

K, N = MWDoubleLayer3D(im*κ), NCross()

e = PlaneWaveMW(z, x, κ, complex(1.0))
H = -1/(im*μ*ω)*curl(E)
h = (n × H) × n

@hilbertspace j; @hilbertspace m;
mfie = @discretise (K+0.5N)[m,j] == h[m]  j∈X m∈Y
u = solve(mfie)

using PlotlyJS
include(Pkg.dir("CompScienceMeshes","examples","matlab_patches.jl"))

fcr, geo = facecurrents(u, X)
p = patch(geo, real.(norm.(fcr)))
