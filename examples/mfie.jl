using CompScienceMeshes
using BEAST
include(Pkg.dir("BEAST","src","lusolver.jl"))

o, x, y, z = euclidianbasis(3)
n = BEAST.n

Γ = readmesh(joinpath(dirname(@__FILE__),"sphere2.in"))
X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

ω = κ = 1.0
μ = ϵ = 1.0

K = MWDoubleLayer3D(im*κ)
N = NCross()

e = PlaneWaveMW(z, x, κ, complex(1.0))
h = -1/(im*μ*ω)*curl(e)
f = (n × h) × n

j, = hilbertspace(:j)
m, = hilbertspace(:m)

MFIE = @varform (K+0.5N)[m,j] == f[m]
mfie = @discretise MFIE j∈X m∈Y

u = solve(mfie)
fcr, geo = facecurrents(u, X)

include(Pkg.dir("CompScienceMeshes","examples","matlab_patches.jl"))
A = real.(norm.(fcr))
p = patch(Γ, A)
