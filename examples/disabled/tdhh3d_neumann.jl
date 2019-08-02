using CompScienceMeshes
using BEAST

G = meshsphere(1.0, 0.25)

c = 1.0
S = BEAST.HH3DDoubleLayerTDBIO(speed_of_light=c)
Id = BEAST.Identity()

# width, delay, scaling = 8.0, 12.0, 1.0
# gaussian = creategaussian(width, delay, scaling)
# e = BEAST.planewave(point(0,0,1), c, gaussian)

# X = lagrangecxd0(G)
X = lagrangecxd0(G)
Y = lagrangec0d1(G)

Δt, Nt = 0.08, 300
T = timebasisc0d1(Δt, Nt)
U = timebasiscxd0(Δt, Nt)

W = Y ⊗ U # test space
V = X ⊗ T # trial space

# nbd = center(chart(G, first(cells(G))))
# refs = refspace(X)
# vals = refs(nbd)

Z = assemble(S, W, V)
W = assemble(Id ⊗ Id, W, V)

Q1 = assemble(Id, Y, X)
Q2 = assemble(Id, X, Y)

X = duallagrangecxd0(G, boundary(G))
Y = lagrangec0d1(G)

Q1 = assemble(Id, X, Y)
Q2 = assemble(Id, Y, X)

using LinearAlgebra
norm(Q1)
norm(Q1-Q2')

Juno.@enter assemble(Id, Y, X)
