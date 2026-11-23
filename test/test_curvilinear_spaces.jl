using Test
using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra

# The spaces built on the barycentric refinement (and the continuous Lagrange space) on meshes
# of curvilinear triangles. A curvilinear mesh of order 1 has the geometry of the flat mesh it is
# made of, so its spaces have the same coefficients; for higher orders the number of functions is
# that of the flat mesh, and the dual constant Lagrange functions keep unit integral on the
# curved refinement.
@testset "spaces on curvilinear meshes" begin

    CSM = CompScienceMeshes
    flat = meshsphere(1.0, 0.6)
    lin = CSM.CurvilinearMesh(vertices(flat), [CSM.HOSimplexGraph{3,1}(c.indices...) for c in cells(flat)])

    shapes(X) = [[(s.cellid, s.refid, s.coeff) for s in fn] for fn in X.fns]
    spaces = (raviartthomas, buffachristiansen, lagrangec0d1, lagrangecxd0, duallagrangec0d1, duallagrangecxd0)

    for f in spaces
        A, B = f(flat), f(lin)
        @test numfunctions(A) == numfunctions(B)
        @test all(length.(A.fns) .== length.(B.fns))
        @test all(((a[1], a[2]) == (b[1], b[2]) && isapprox(a[3], b[3]; rtol=1e-12))
                  for (fa, fb) in zip(shapes(A), shapes(B)) for (a, b) in zip(fa, fb))
        @test all(isapprox.(A.pos, B.pos; atol=1e-12))
    end

    for O in (2, 3)
        Γ = gmshsphere(1.0, 0.6; order=O)
        @test Γ isa CSM.CurvilinearMesh
        for f in spaces
            @test numfunctions(f(Γ)) == numfunctions(f(flat))
        end

        # the Buffa-Christiansen functions are supported on the curved refinement, and positioned
        # at the (curved) midpoints of the edges, which lie on the sphere
        Y = buffachristiansen(Γ)
        @test Y.geo isa CSM.BarycentricRefinement
        @test Y.geo.mesh isa CSM.CurvilinearMesh
        @test all(p -> isapprox(norm(p), 1.0; atol=1e-2), Y.pos)

        X = duallagrangecxd0(Γ)
        fine = X.geo
        @test all(fn -> isapprox(sum(s.coeff * volume(chart(fine, s.cellid)) for s in fn), 1; atol=1e-12), X.fns)
    end
end
