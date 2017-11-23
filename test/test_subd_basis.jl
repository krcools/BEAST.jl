using CompScienceMeshes,BEAST
# fn = joinpath(dirname(@__FILE__),"sphere.in")
# G = readmesh(fn)
G0 = meshsphere(1.0, 0.5)
G = Loop_subdivision(G0)
X = subdsurface(G)

els, ad = BEAST.assemblydata(X)

identityop    = Identity()
singlelayer = BEAST.HH3DSingleLayerFDBIO(1.0)
I = assemble(identityop, X, X)
S = assemble(singlelayer, X, X)
ncd = cond(I)
