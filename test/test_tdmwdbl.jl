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
δ = timebasisdelta(Δt, Nt)
T = timebasiscxd0(Δt, Nt)
T2 = timebasisshiftedlagrange(Δt, Nt, 2)
V = X ⊗ T
W = Y ⊗ T

Z, store1 = BEAST.allocatestorage(K, W, V,
    Val{:densestorage}, BEAST.LongDelays{:ignore})

BEAST.assemble!(K, W, V, store1)
@test all(==(0), Z[:,:,1])

W = X⊗δ
V = Y⊗T2

K = TDMaxwell3D.doublelayer(speedoflight=1.0, numdiffs=1)
Z2, store2 = BEAST.allocatestorage(K, W, V, Val{:densestorage}, BEAST.LongDelays{:ignore})
BEAST.assemble!(K, W, V, store2)
@test all(==(0), Z2[:,:,1])

γ = geometry(Y)
ch1 = chart(G, first(cells(G)))
qps = quadpoints(ch1, 3)

x = qps[1][1]
w = qps[1][2]

test_els = BEAST.elements(G)
trial_els = BEAST.elements(γ)
time_els = [(0.0, 20.0)]

qdt = BEAST.quaddata(K, refspace(X), refspace(Y), refspace(T2), test_els, trial_els, nothing)
qrl = BEAST.quadrule(K,
    refspace(X), refspace(Y), refspace(T2),
    1, test_els[1], 1, trial_els[1], 1, time_els[1], qdt)
z = zeros(Float64, 3, 3, 10)
BEAST.innerintegrals!(z, K, x, refspace(X), refspace(Y), refspace(T2), test_els[1], trial_els[1], (0.0, 20.0), qrl, w)
@test_broken !any(isnan.(z))

verts = trial_els[1].vertices
import WiltonInts84
iG, ivG, igG = WiltonInts84.wiltonints(verts[1], verts[2], verts[3], cartesian(x), 0.0, 20.0, Val{2}, qrl.workspace)
@test_broken !any(isnan(iG))

ctr = WiltonInts84.contour!(verts[1],verts[2],verts[3],cartesian(x),0.0,20.0,qrl.workspace)
iG, ivG, igG = WiltonInts84.wiltonints(ctr,cartesian(x),Val{2})
@test_broken !any(isnan(iG))

import WiltonInts84
a = -0.23570226039551564
b = 1.5700924586837734e-16
p = -6.371760312437358e-16
h = 0.0
using StaticArrays
m = SVector(-0.7071067811865474, 0.7071067811865478, -0.0)
iG, vG = WiltonInts84.segintsg(a, b, p, h, m, Val{2})
@test_broken !any(isnan.(iG))
