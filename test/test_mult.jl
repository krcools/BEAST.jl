using CompScienceMeshes
using BEAST

using Test
using LinearAlgebra

for T in [Float32, Float64]
faces = meshrectangle(T(1.0), T(1.0), T(0.5), 3)
srt_bnd_faces = sort.(boundary(faces))
local edges = submesh(skeleton(faces,1)) do edge
    !(sort(edge) in srt_bnd_faces)
end

srt_bnd_nodes = sort.(skeleton(boundary(faces),0))
@test length(srt_bnd_nodes) == 8
nodes = submesh(skeleton(faces,0)) do node
    !(sort(node) in srt_bnd_nodes)
end
@test length(nodes) == 1

Conn = connectivity(nodes, edges, sign)

local X = raviartthomas(faces, cellpairs(faces,edges))
@test numfunctions(X) == 8

divX = divergence(X)
Id = BEAST.Identity()
DD = assemble(Id, divX, divX)
@test rank(DD) == 7
L = divX * Conn
for sh in L.fns[1]
    @test isapprox(sh.coeff, 0, atol=1e-8)
end
end