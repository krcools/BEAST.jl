using Test

using BEAST
using CompScienceMeshes

o, x, y = point(0,0,0), point(1,0,0), point(0,1,0)
G = Mesh([o,x,y,x+y],[CompScienceMeshes.SimplexGraph(1,2,3)])

@test numvertices(G) == 4
@test numcells(G) == 1

c = 1.0
K = BEAST.HH3DDoubleLayerTDBIO(speed_of_light=c)

X = lagrangecxd0(G)
@test numfunctions(refspace(X), domain(chart(G, first(G)))) == 1

Y = duallagrangec0d1(G)
@test numfunctions(Y) == 1

Gref = geometry(Y)
@test numfunctions(refspace(Y), domain(chart(Gref, first(Gref)))) == 3

# X = duallagrangecxd0(G, boundary(G))
# Y = lagrangec0d1(G, dirichlet=false)
# @test numfunctions(X) == 3
# @test numfunctions(Y) == 3

Δt = 20.0
Nt = 5
T = timebasiscxd0(Δt, Nt)
V = X ⊗ T
W = Y ⊗ T

fr2, store2 = BEAST.allocatestorage(K, V, V, Val{:densestorage}, BEAST.LongDelays{:ignore})
BEAST.assemble!(K, W, V, store2)
Z = fr2()

T = scalartype(K,V,V)
Z1 = zeros(T, size(Z)[1:2])
BEAST.ConvolutionOperators.timeslice!(Z1, Z, 1)
@test all(==(0), Z1)

# import WiltonInts84
# trial_element = chart(G, first(cells(G)))
# x = neighborhood(trial_element,(1/3,1/3))
# x = neighborhood(trial_element,(0.0,0.0))
# r, R = 0.0, 20.0
# workspace = WiltonInts84.workspace(eltype(τ.vertices))
# G, vG, GG = WiltonInts84.wiltonints(
#     trial_element[1],
#     trial_element[2],
#     trial_element[3],
#     cartesian(x), r, R, Val{0},workspace)
