using CompScienceMeshes
using BEAST

sphere = CompScienceMeshes.tetmeshsphere(1.0, 0.35)
hemi = submesh(tet -> cartesian(center(chart(sphere,tet)))[3] < 0, sphere)

X = BEAST.nedelecd3d(sphere)
Y = BEAST.restrict(X, hemi)
