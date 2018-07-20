using CompScienceMeshes
using BEAST
using Test
G = readmesh(joinpath(dirname(@__FILE__),"assets","sphere872.in"))
 # G0 = readmesh(fn)
# G0 = meshsphere(1.0, 0.5)
#G = readmesh("/Users/Benjamin/Documents/sphere.in")
# G = Loop_subdivision(G0)
X = subdsurface(G)

#els, ad = BEAST.assemblydata(X)

identityop  = Identity()
singlelayer = Helmholtz3D.singlelayer(gamma=1.0)
I = assemble(identityop, X, X)
#S = assemble(singlelayer, X, X)
ncd = cond(I)

@test ncd ≈ 64.50401358713235 rtol=1e-8
