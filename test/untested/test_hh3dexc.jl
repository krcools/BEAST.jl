using CompScienceMeshes
using BEAST
using Test

sphere = readmesh(joinpath(dirname(@__FILE__),"assets","sphere3.in"))

κ = 2π
direction = point(0,0,1)
f = BEAST.HH3DPlaneWave(direction, κ)

v1 = f(point(0,0,0))
v2 = f(point(0,0,0.5))

@test v1 ≈ +1
@test v2 ≈ -1

import BEAST.∂n
p = ∂n(f)

s = chart(m,1)
c = neighborhood(s, [1,1]/3)

r = cartesian(c)
n = normal(s)

w1 = p(c)
w2 = -im*κ*dot(direction, n)*f(r)

w1 ≈ w2

N = BEAST.HH3DHyperSingularFDBIO(im*κ)
X = BEAST.lagrangec0d1(m)

numfunctions(X)

Nxx = assemble(N, X, X)
