using CompScienceMeshes,BEAST
G = readmesh("/Users/Benjamin/Documents/sphere.in")
X = subdsurface(G)
identityop    = Identity()
singlelayer = BEAST.HH3DSingleLayerFDBIO(1.0)
I = assemble(identityop, X, X)
S = assemble(singlelayer, X, X)
ncd = cond(I)
