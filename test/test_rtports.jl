using Test
using BEAST
using CompScienceMeshes

meshZ = meshrectangle(1.0,1.0,0.5)
@test numvertices(meshZ) == 9

translate!(meshZ,point(0.0,0.0,0.5))
mesh0 = meshrectangle(1.0,1.0,0.5)
@test numvertices(mesh0) == 9

idcsZ = meshZ.faces[1]
@test size(idcsZ.indices) == (3,)

# vertsZ = vertices(meshZ, idcsZ)
vertsZ = chart(meshZ, idcsZ).vertices
@test size(vertsZ) == (3,)

Γ = meshsegment(1.0,0.5,3)
@test numvertices(Γ) == 3

translate!(Γ,point(0.0,0.0,0.5))
γ = meshsegment(1.0,0.5,3)

@test numvertices(γ) == 3

idcsΓ = Γ.faces[1]
@test size(idcsΓ.indices) == (2,)

# vertsΓ = vertices(Γ, idcsΓ)
vertsΓ = chart(Γ, idcsΓ).vertices
@test size(vertsΓ) == (2,)

sum_mesh = weld(meshZ,mesh0)
@test numvertices(sum_mesh) == 18

edges = skeleton(sum_mesh,1)
@test numcells(edges) == 32

cps = cellpairs(sum_mesh, edges)
@test size(cps) == (2,32)

# select only outer edges
overlaps = overlap_gpredicate(Γ)
on_junction = (m,c) -> overlaps(chart(m,c))

# pred = x -> on_junction(x)
edges = submesh(on_junction, skeleton(meshZ,1))
# edges = skeleton(pred, meshZ, 1)
@test length(edges) == 2
# @test size(edges.faces[1]) == (2,)
@test dimension(edges) == 1

cps = cellpairs(sum_mesh, edges)
@test size(cps) == (2,2)

# test portcells function
pc = portcells(sum_mesh, Γ)
@test size(pc) == (2,2)

# #test rt_vedge function
# ce1 = rt_vedge(pc,1.0)
# @test size(ce1) == (1,)
#
# #test rt_cedge function
# ce1 = rt_cedge(pc,1.0)
# @test size(ce1) == (2,)

# build the Raviart-Thomas elements
rt = rt_ports(sum_mesh, [Γ, γ])
@test numfunctions(rt) == 19

rt = rt_ports(sum_mesh, [Γ, γ], [Γ, γ]) #Don't do this normally.
@test numfunctions(rt) == 22      #Only used here to check numfunctions

#
