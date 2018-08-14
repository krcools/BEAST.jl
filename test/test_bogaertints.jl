using Test
using CompScienceMeshes
using StaticArrays

import BEAST
BE = BEAST

Γ = readmesh(joinpath(dirname(@__FILE__),"assets","sphere2.in"))
nc = numcells(Γ)
t = chart(Γ, Γ.faces[1])
s = chart(Γ, Γ.faces[nc])

X = BEAST.raviartthomas(Γ)
x = BEAST.refspace(X)

κ = 1.0
op = BEAST.MWSingleLayer3D(κ)

n = BE.numfunctions(x)
z1 = zeros(ComplexF64, n, n)
z2 = zeros(z1)

tqd = BE.quadpoints(x, [t], (12,13))
bqd = BE.quadpoints(x, [t], (13,))

SE_strategy = BE.WiltonSEStrategy(
  tqd[2,1],
  BE.DoubleQuadStrategy(
	tqd[1,1],
	bqd[1,1],
  ),
)
BEAST.momintegrals!(op, x, x, t, t, z1, SE_strategy)

EE_strategy = BEAST.BogaertSelfPatchStrategy(20)
BEAST.momintegrals!(op, x, x, t, t, z2, EE_strategy)

Γ = meshrectangle(1.0, 1.0, 1.0)
t = chart(Γ, Γ.faces[1])
s = chart(Γ, Γ.faces[2])


z3 = zeros(ComplexF64, n, n)
EE_strategy = BEAST.BogaertEdgePatchStrategy(13,30)
BEAST.momintegrals!(op, x, x, t, s, z3, EE_strategy)

z4 = zeros(ComplexF64, n, n)
tqd = BE.quadpoints(x, [t], (12,13))
bqd = BE.quadpoints(x, [s], (13,))
SE_strategy = BE.WiltonSEStrategy(
  tqd[2,1],
  BE.DoubleQuadStrategy(
	tqd[1,1],
	bqd[1,1],
  ),
)
BEAST.momintegrals!(op, x, x, t, s, z4, SE_strategy)

p = [
    point(0.0, 0.0, 0.0),
    point(0.0, 1.0, 0.0),
    point(1.0, 1.0, 0.0),
    point(0.0, -1.0, 0.0),
    point(1.0, -1.0, 0.0)]

t = simplex(p[1], p[2], p[3])
s = simplex(p[1], p[4], p[5])

z5 = zeros(ComplexF64, n, n)
EE_strategy = BEAST.BogaertPointPatchStrategy(9,10)
BEAST.momintegrals!(op, x, x, t, s, z5, EE_strategy)

z6 = zeros(ComplexF64, n, n)
tqd = BE.quadpoints(x, [t], (12,13))
bqd = BE.quadpoints(x, [s], (13,))
SE_strategy = BE.WiltonSEStrategy(
  tqd[2,1],
  BE.DoubleQuadStrategy(
	tqd[1,1],
	bqd[1,1],
  ),
)
BEAST.momintegrals!(op, x, x, t, s, z6, SE_strategy)

@test norm(z1-z2)/norm(z1) < 1.0e-6
@test norm(z3-z4)/norm(z3) < 3.0e-6
@test norm(z5-z6)/norm(z5) < 1.0e-8

@show norm(z1-z2)/norm(z1)
@show norm(z3-z4)/norm(z3)
@show norm(z5-z6)/norm(z5)

# repeat the test for the MFIE interactions
κ = 1.0
op = BEAST.MWDoubleLayer3D(0.0)

p = [
    point(0.0, 0.0, 0.0),
    point(1.0, 0.0, 0.0),
    point(1.0, 1.0, 0.0),
    point(0.0, 1.0, 1.0),
    point(1.0,-1.0, 0.0)
]

t1 = simplex(p[1], p[2], p[3])
t2 = simplex(p[1], p[3], p[4])
t3 = simplex(p[1], p[2], p[5])


# test the edge patch integral
z1 = zeros(ComplexF64, BE.numfunctions(x), BE.numfunctions(x))
z2 = zeros(ComplexF64, BE.numfunctions(x), BE.numfunctions(x))
z3 = zeros(ComplexF64, BE.numfunctions(x), BE.numfunctions(x))

s1 = BEAST.BogaertEdgePatchStrategy(13, 30)
tqd = BE.quadpoints(x, [t1], (12,13))
bqd = BE.quadpoints(x, [t2], (13,))
s2 = BE.WiltonSEStrategy(
  tqd[2,1],
  BE.DoubleQuadStrategy(
	tqd[1,1],
	bqd[1,1],
  ),
)
s3 = BE.DoubleQuadStrategy(
  tqd[1,1],
  tqd[1,1],
)

BEAST.momintegrals!(op, x, x, t1, t2, z1, s1)
BEAST.momintegrals!(op, x, x, t1, t2, z2, s2)
BEAST.momintegrals!(op, x, x, t1, t2, z3, s3)

@test norm(z2-z1)/norm(z2) < 1.0e-3

s1 = BEAST.BogaertPointPatchStrategy(9,10)
tqd = BE.quadpoints(x, [t2], (12,13))
bqd = BE.quadpoints(x, [t3], (13,))
s2 = BE.WiltonSEStrategy(
  tqd[2,1],
  BE.DoubleQuadStrategy(
	tqd[1,1],
	bqd[1,1],
  ),
)
s3 = BE.DoubleQuadStrategy(
  tqd[1,1],
  tqd[1,1],
)

fill!(z1, 0); BEAST.momintegrals!(op, x, x, t2, t3, z1, s1)
fill!(z2, 0); BEAST.momintegrals!(op, x, x, t2, t3, z2, s2)
fill!(z3, 0); BEAST.momintegrals!(op, x, x, t2, t3, z3, s3)

@test norm(z2-z1)/norm(z2) < 2.0e-6

s1 = BEAST.BogaertSelfPatchStrategy(20)
tqd = BE.quadpoints(x, [t1], (12,13))
bqd = BE.quadpoints(x, [t1], (13,))
s2 = BE.WiltonSEStrategy(
  tqd[2,1],
  BE.DoubleQuadStrategy(
	tqd[1,1],
	bqd[1,1],
  ),
)
s3 = BE.DoubleQuadStrategy(
  tqd[1,1],
  tqd[1,1],
)

fill!(z1, 0); BEAST.momintegrals!(op, x, x, t1, t1, z1, s1)
fill!(z2, 0); BEAST.momintegrals!(op, x, x, t1, t1, z2, s2)
fill!(z3, 0); BEAST.momintegrals!(op, x, x, t1, t1, z3, s3)



# Test if IBPP gives the correct value for ∫∫∇G
t2 = simplex(p[3], p[4], p[1])
t3 = simplex(p[2], p[5], p[1])

s1 = BEAST.BogaertPointPatchStrategy(9, 10)
_, I = BEAST.GetIntegrals(t2[1],t2[2],t2[3],t3[1],t3[2], 0.0, s1)
#I = I.GG

Rpnt = eltype(p)
Cplx = typeof(complex(zero(eltype(p[1]))))
Cpnt = similar_type(Rpnt, Cplx, Size(length(Rpnt),))

u, w = BEAST.trgauss(10)
J = zeros(Cpnt, 3, 3)
for i in 1:length(w)
    p = neighborhood(t2, u[:,i])
    x = cartesian(p)
    λ = barycentric(p)
    dx = w[i] * jacobian(p)
    for j in 1:length(w)
        q = neighborhood(t3, u[:,j])
        y = cartesian(q)
        μ = barycentric(q)
        dy = w[j] * jacobian(q)
        r = norm(x-y)
        for i in 1:3
            for j in 1:3
                J[i,j] += dx * dy / r^3 * (x-y) / 4π * λ[i] * μ[j]
            end
        end
    end
end

@test maximum([real(norm(x-y)) for (x,y) in zip(I,J)]) < 2.0e-6
