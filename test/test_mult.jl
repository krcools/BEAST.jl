using CompScienceMeshes
using BEAST

using Test
using LinearAlgebra

for T in [Float32, Float64]
    faces = meshrectangle(T(1.0), T(1.0), T(0.5), 3)

    local bnd = boundary(faces)
    local edges = submesh(!in(bnd), skeleton(faces,1))

    local bnd_nodes = skeleton(bnd, 0)
    @test length(bnd_nodes) == 8
    local nodes = submesh(!in(bnd_nodes), skeleton(faces,0))
    @test length(nodes) == 1

    Conn = connectivity(nodes, edges, sign)

    local X = raviartthomas(faces, cellpairs(faces,edges))
    @test numfunctions(X) == 8

    divX = divergence(X)
    Id = BEAST.Identity()
    DD = Matrix(assemble(Id, divX, divX))
    @test rank(DD) == 7
    L = divX * Conn
    for sh in L.fns[1]
        @test isapprox(sh.coeff, 0, atol=1e-8)
    end
end