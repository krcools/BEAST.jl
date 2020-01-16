using CompScienceMeshes
using BEAST

m = CompScienceMeshes.tetmeshsphere(1.0, 0.45)
bnd_m = boundary(m)

m1 = skeleton(m,1)
X = BEAST.nedelecc3d(m,m1)

Y = BEAST.ttrace(X,bnd_m)
@assert length(geometry(Y)) == length(bnd_m)
