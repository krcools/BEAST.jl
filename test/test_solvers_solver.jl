using Test
using BEAST
using CompScienceMeshes
using LinearAlgebra

fn = joinpath(pkgdir(BEAST), "test", "assets", "sphere45.in")
m = CompScienceMeshes.readmesh(fn)

X = raviartthomas(m)
𝕏 = X × X

I = BEAST.Identity()
T = Maxwell3D.singlelayer(wavenumber=1.0)

@hilbertspace m j
@hilbertspace k l

a = (
    I[k,m] + 2*T[k,j] +
    3im*T[l,m] - I[l,j]
)

A = assemble(a, 𝕏, 𝕏; threading=:cellcoloring)
# this test is brittle but not sure how to do this otherwise...
@test A.maps[1].lmap.A.lmap === A.maps[4].lmap.A.lmap 
@test A.maps[2].lmap.A.lmap === A.maps[3].lmap.A.lmap