@testitem "Lagrange spaces do not depend on the unit of length" begin
    using CompScienceMeshes
    using LinearAlgebra

    # The global degrees of freedom identify the vertices of a node or an edge with those
    # of a cell by their positions. With an absolute tolerance of sqrt(eps), all
    # vertices of a cell smaller than about 1e-8 were identified with its first one,
    # and the construction failed with a singular matrix.
    scaled(mesh, s) = Mesh([s * p for p in vertices(mesh)], copy(cells(mesh)))
    shapes(space, i) = sort([(f.cellid, f.refid, f.coeff) for f in space.fns[i]])
    same_shapes(a, b) =
        length(a) == length(b) &&
        all(p[1:2] == q[1:2] && isapprox(p[3], q[3]; rtol=1e-12) for (p, q) in zip(a, b))

    # The edges of skeleton(mesh, 1) are sorted along a space-filling curve of their
    # coordinates, so the numbering of the edge functions (order 2) can change with the
    # scale. The functions are therefore matched by their positions.
    function matching(space, reference, s)
        return map(eachindex(reference.pos)) do i
            j = findfirst(p -> isapprox(p, s * reference.pos[i]; rtol=1e-12), space.pos)
            return something(j, 0)
        end
    end

    mesh = meshrectangle(2.0, 1.0, 0.25)
    for order in (1, 2), dirichlet in (false, true)
        reference = lagrangec0(mesh; dirichlet, order)
        G = assemble(Identity(), reference, reference)
        for s in (1e-12, 1e-9, 1e-6, 1e3)
            space = lagrangec0(scaled(mesh, s); dirichlet, order)
            @test numfunctions(space) == numfunctions(reference)
            p = matching(space, reference, s)
            @test sort(p) == 1:numfunctions(reference)
            # The same functions: the same shapes on the same cells, and a Gram matrix
            # that scales with the area
            @test all(same_shapes(shapes(space, p[i]), shapes(reference, i)) for i in eachindex(p))
            @test assemble(Identity(), space, space)[p, p] ≈ s^2 * G rtol = 1e-10
        end
    end

    # The space of lagrangec0d1, which the other bases build on
    for dirichlet in (false, true)
        reference = lagrangec0d1(mesh; dirichlet)
        space = lagrangec0d1(scaled(mesh, 1e-9); dirichlet)
        @test all(same_shapes(shapes(space, i), shapes(reference, i)) for i in eachindex(space.fns))
    end
end
