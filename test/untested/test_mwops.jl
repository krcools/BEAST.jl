using CompScienceMeshes
using BEAST
using Test

############################################################################################
#
# Test the sign of the double layer potential
#
############################################################################################

γ = 0.0
K = MWDoubleLayer3D(γ)

τ = Mesh(
    [
        point(-4,0,0),
        point(-6,0,0),
        point(-5,-1,0)
    ],
    [
        index(1,2,3)])

σ = Mesh(
    [
        point(4,0,0),
        point(6,0,0),
        point(5,0,-1)
    ],
    [
        index(1,2,3)
    ]
)

@test numcells(τ) == 1
@test numcells(σ) == 1

T = raviartthomas(τ, boundary(τ))
S = raviartthomas(σ, boundary(σ))

@test numfunctions(T) == 3
@test numfunctions(S) == 3

t = chart(τ, first(cells(τ)))
s = chart(σ, first(cells(σ)))

p = neighborhood(t, [1,1]/2)
q = neighborhood(s, [1,1]/2)

R = norm(cartesian(p)-cartesian(q))

rt = refspace(T)
f = (refspace(T))(p)[3][1]
g = (refspace(S))(q)[3][1]

@test dot(f, point(0,1,0)) > 0
@test dot(g, point(0,0,1)) > 0

@show 1/R^2

z = assemble(K, T, S)

qd = quaddata(K, rt, rt, elements(τ), elements(σ))
zlocal = zeros(promote_type(scalartype(K), scalartype(S), scalartype(T)), 3, 3)
strat = BEAST.quadrule(K, rt, rt, 1, t, 1, s, qd)
strat =  BEAST.DoubleQuadStrategy(qd.tpoints[1,1],  qd.bpoints[1,1])
BEAST.momintegrals!(K, rt, rt, t, s, zlocal, strat)
z = zlocal[3,3]
