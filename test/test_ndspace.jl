using CompScienceMeshes
using BEAST
using Test

fn = joinpath(dirname(@__FILE__),"assets","rect1.in")

mesh = readmesh(fn)
edgs = skeleton(mesh,1)

charts = [chart(edgs,edg) for edg in cells(edgs)]
ctrs = [cartesian(center(cht)) for cht in charts]

ND = BEAST.nedelec(mesh, edgs)

@test numfunctions(ND) == numcells(edgs)
@test length(ND.fns[1]) == 1
@test length(ND.fns[2]) == 1
@test length(ND.fns[3]) == 2
@test length(ND.fns[4]) == 1
@test length(ND.fns[5]) == 1

# center of the diagonal
ctr = center(charts[3])
face_charts = [chart(mesh,fce) for fce in cells(mesh)]
nbd1 = neighborhood(face_charts[1], carttobary(face_charts[1], ctr))
nbd2 = neighborhood(face_charts[2], carttobary(face_charts[2], ctr))
t = tangents(ctr,1)
ut = t / norm(t)

ndlocal = refspace(ND)
fn = ND.fns[3]
v1 = dot(fn[1].coeff*ndlocal(nbd1)[1][1],ut)
v2 = dot(fn[2].coeff*ndlocal(nbd2)[2][1],ut)
@test v1 ≈ v2 ≈ -1/√2 ≈ -1/volume(charts[3])
