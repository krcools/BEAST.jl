using Test
using BEAST
using CompScienceMeshes
using LinearAlgebra

# Local operators between the spaces on a curvilinear mesh and those on its barycentric
# refinement, in both directions: the test geometry refining the trial geometry and conversely.
# Identity is symmetric and NCross antisymmetric, so the two directions are each other's
# transposes. A curvilinear mesh of order 1 has the geometry of the flat mesh it is made of, so
# its matrices are those of the flat mesh.
@testset "local operators on curvilinear refinements" begin

    CSM = CompScienceMeshes
    asm(op, X, Y) = Matrix(assemble(op, X, Y; verbose=false))
    pairs(Γ) = ((raviartthomas(Γ), buffachristiansen(Γ)), (lagrangecxd0(Γ), duallagrangec0d1(Γ)))
    ops = ((BEAST.Identity(), 1), (BEAST.NCross(), -1))

    flat = meshsphere(1.0, 0.6)
    lin = CSM.CurvilinearMesh(vertices(flat), [CSM.HOSimplexGraph{3,1}(c.indices...) for c in cells(flat)])

    for ((Xf, Yf), (Xl, Yl)) in zip(pairs(flat), pairs(lin)), (op, _) in ops
        op isa BEAST.NCross && !(Xf isa BEAST.RTBasis) && continue
        A, B = asm(op, Xl, Yl), asm(op, Xf, Yf)
        @test norm(A - B) ≤ 1e-12 * norm(B)
    end

    for O in (2, 3)
        Γ = gmshsphere(1.0, 0.6; order=O)
        for (X, Y) in pairs(Γ), (op, s) in ops
            op isa BEAST.NCross && !(X isa BEAST.RTBasis) && continue
            A, B = asm(op, X, Y), asm(op, Y, X)
            @test norm(A - s * transpose(B)) ≤ 1e-12 * norm(A)
        end
    end
end
