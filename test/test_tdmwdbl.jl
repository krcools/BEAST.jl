using Test

using BEAST
using CompScienceMeshes

o, x, y = point(0,0,0), point(1,0,0), point(0,1,0)
G = Mesh([o,x,y],[index(1,2,3)])

@test numvertices(G) == 3
@test numcells(G) == 1

K = MWDoubleLayerTDIO(1.0, 1.0, 0)

X = raviartthomas(G, skeleton(G,1))
Y = buffachristiansen(G, boundary(G))
@test numfunctions(X) == 3
@test numfunctions(Y) == 3

Δt = 20.0
Nt = 5
T = timebasiscxd0(Δt, Nt)
V = X ⊗ T
W = Y ⊗ T

Z, store1 = BEAST.allocatestorage(K, W, V, Val{:densestorage})
BEAST.assemble!(K, W, V, store1)
@test all(==(0), Z[:,:,1])
