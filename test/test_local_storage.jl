using Test
using BEAST
using CompScienceMeshes

fn = joinpath(@__DIR__, "assets/sphere5.in")
m = readmesh(fn)

Id = BEAST.Identity()
X = BEAST.raviartthomas(m)

Z1 = assemble(Id, X, X, storage_policy=Val{:densestorage})
Z2 = assemble(Id, X, X, storage_policy=Val{:bandedstorage})

@test Z1 ≈ Z2 atol=1e-8