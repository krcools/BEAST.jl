using CompScienceMeshes
using BEAST
using LinearAlgebra

G = meshsphere(1.0, 0.25)

c = 1.0
S = BEAST.HH3DSingleLayerTDBIO(c)
D = BEAST.HH3DDoubleLayerTDBIO(speed_of_light=c)
Id = BEAST.Identity()

width, delay, scaling = 8.0, 12.0, 1.0
gaussian = creategaussian(width, delay, scaling)
e = BEAST.planewave(point(0,0,1), c, gaussian)
de = BEAST.planewave(point(0,0,1), c, derive(gaussian))
h = BEAST.gradient(e)

# X = lagrangecxd0(G)
X = lagrangecxd0(G)
Y = duallagrangec0d1(G)
# Y = lagrangecxd0(G)

# X = duallagrangecxd0(G, boundary(G))
# Y = lagrangec0d1(G)

Δt, Nt = 0.08, 300
# Δt, Nt = 0.16, 300
# T = timebasisc0d1(Δt, Nt)
P = timebasiscxd0(Δt, Nt)
H = timebasisc0d1(Δt, Nt)
δ = timebasisdelta(Δt, Nt)

# assemble the right hand side

bd = assemble(n⋅h, Y ⊗ P)
Z1d = assemble(Id ⊗ Id, Y ⊗ P, X ⊗P)
Z0d = assemble(D, Y ⊗ P, X ⊗P)
Zd = Z0d + (0.5)*Z1d
u = marchonintime(inv(Z[:,:,1]), Zd, bd, Nt)

bs = assemble(e, X ⊗ δ)
Zs = assemble(S, X ⊗ δ, X ⊗ H)
v = marchonintime(inv(Zs[:,:,1]), Zs, bs, Nt)
