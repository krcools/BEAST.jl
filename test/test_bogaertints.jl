using FactCheck
using CompScienceMeshes
using StaticArrays

import BEAST
BE = BEAST

Γ = meshsphere(1.0, 0.2)
nc = numcells(Γ)
#t = simplex(cellvertices(Γ,1))
t = simplex(vertices(Γ, Γ.faces[1]))
#s = simplex(cellvertices(Γ,nc))
s = simplex(vertices(Γ, Γ.faces[nc]))

X = BEAST.raviartthomas(Γ)
x = BEAST.refspace(X)

κ = 1.0
op = BEAST.MWSingleLayer3D(κ)

n = BE.numfunctions(x)
z1 = zeros(Complex128, n, n)
z2 = zeros(z1)

tqd = BE.quadpoints(x, [t], (12,13))
bqd = BE.quadpoints(x, [t], (13,))

#SE_strategy = BEAST.WiltonSEStrategy(13,12,13)
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
#t = simplex(cellvertices(Γ,1))
t = simplex(vertices(Γ, Γ.faces[1]))
#s = simplex(cellvertices(Γ,2))
s = simplex(vertices(Γ, Γ.faces[2]))


z3 = zeros(Complex128, n, n)
EE_strategy = BEAST.BogaertEdgePatchStrategy(13,30)
BEAST.momintegrals!(op, x, x, t, s, z3, EE_strategy)

z4 = zeros(Complex128, n, n)
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

z5 = zeros(Complex128, n, n)
EE_strategy = BEAST.BogaertPointPatchStrategy(9,10)
BEAST.momintegrals!(op, x, x, t, s, z5, EE_strategy)

z6 = zeros(Complex128, n, n)
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

facts("Test IB Self/Edge/Point-patch strategy...") do
    @fact norm(z1-z2)/norm(z1) --> less_than(1.0e-6)
    @fact norm(z3-z4)/norm(z3) --> less_than(3.0e-6)
    @fact norm(z5-z6)/norm(z5) --> less_than(1.0e-8)
end

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
z1 = zeros(Complex128, BE.numfunctions(x), BE.numfunctions(x))
z2 = zeros(Complex128, BE.numfunctions(x), BE.numfunctions(x))
z3 = zeros(Complex128, BE.numfunctions(x), BE.numfunctions(x))

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

facts("Test IB Edge-patch for MFIE interactions...") do
    @fact norm(z2-z1)/norm(z2) --> less_than(1.0e-3)
end

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

facts("Test IB Point-patch for MFIE interactions...") do
    @fact norm(z2-z1)/norm(z2) --> less_than(2.0e-6)
end

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

facts("Test IB Point-patch primitve integrals") do
        @fact maximum([real(norm(x-y)) for (x,y) in zip(I,J)]) --> less_than(2.0e-6)
end
