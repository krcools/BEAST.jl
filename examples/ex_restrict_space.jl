using CompScienceMeshes
using BEAST

sphere = CompScienceMeshes.tetmeshsphere(1.0, 0.35)
hemi = submesh((m,tet) -> cartesian(CompScienceMeshes.center(chart(m,tet)))[3] < 0, sphere)

X = BEAST.nedelecd3d(sphere)
Y = BEAST.restrict(X, hemi)
