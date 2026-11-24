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

@testitem "Lagrange spaces distinguish the vertices of thin cells" begin
    using CompScienceMeshes
    using StaticArrays

    function triangle(T, scale, aspect)
        points = [
            SVector{3,T}(0, 0, 0),
            SVector{3,T}(scale, 0, 0),
            SVector{3,T}(0, scale * aspect, 0),
        ]
        return Mesh(points, [CompScienceMeshes.SimplexGraph(1, 2, 3)])
    end
    # On one triangle each nodal or edge function has one local shape. Their local
    # coefficients are unchanged by an invertible affine map of the triangle; the
    # global edge order may change, so compare by local reference index.
    signature(basis) = sort([(only(fn).refid, Float64(only(fn).coeff)) for fn in basis.fns])

    for (T, aspect, scales) in (
        (Float64, 1e-9, (1e-12, 1e-6, 1.0, 1e3)),
        (Float32, 1e-4, (1e-6, 1.0, 1e3)),
        (Float16, 0.02, (0.1, 1.0, 10.0)),
    ), order in (1, 2)
        reference = lagrangec0(triangle(T, 1.0, 1.0); order, dirichlet=false)
        for scale in scales
            @testset "$T, order=$order, scale=$scale" begin
                basis = lagrangec0(triangle(T, scale, aspect); order, dirichlet=false)
                @test numfunctions(basis) == numfunctions(reference)
                @test all(length(fn) == 1 for fn in basis.fns)
                @test first.(signature(basis)) == first.(signature(reference))
                @test last.(signature(basis)) ≈ last.(signature(reference)) rtol=8eps(T)
                @test all(p -> all(isfinite, p), basis.pos)
            end
        end
    end
end

@testitem "Trace vertex matching is geometric and one-to-one" begin
    using CompScienceMeshes
    using StaticArrays

    # Construct the edge independently from the face, with a small perturbation of
    # its coordinates. For these thin cells the tolerance also includes the wrong
    # vertex: choosing the first acceptable vertex instead of the nearest collapses
    # the edge. The expected pullback follows the intended face edge and its order.
    for (T, aspect) in ((Float64, 1e-9), (Float32, 1e-4), (Float16, 0.02)),
        scale in (0.1, 1.0, 10.0)

        points = [
            SVector{3,T}(0, 0, 0),
            SVector{3,T}(scale, 0, 0),
            SVector{3,T}(0, scale * aspect, 0),
        ]
        face = simplex(points...)
        delta = T(scale * aspect / 100) * SVector{3,T}(1, -1, 0)
        for indices in ((3, 2), (2, 3))
            a, b = points[indices[1]], points[indices[2]]
            expected = simplex(a, b)
            independent = simplex(a + delta, b - delta)
            traced = BEAST.trace(independent, face, cartesian)
            for s in T.((0, 0.3, 1))
                @test traced(s) ≈ cartesian(neighborhood(expected, s)) rtol=8eps(T)
            end
        end
    end

    @testset "Invalid correspondences fail explicitly" begin
        points = [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(0.0, 0.02, 0.0)]
        match(vs) = BEAST.trace_vertex_permutation(vs, points, 0.03)
        @test collect(match(copy(points))) == [1, 2, 3]
        # Equidistant from two vertices, outside the tolerance, or two edge vertices
        # mapping to one face vertex; duplicate face vertices have no unique identity.
        @test_throws ArgumentError match([SVector(0.0, 0.01, 0.0)])
        @test_throws ArgumentError match([SVector(2.0, 0.0, 0.0)])
        @test_throws ArgumentError match([points[1], SVector(0.0, 1e-6, 0.0)])
        @test_throws ArgumentError BEAST.trace_vertex_permutation(
            [points[1]], [points[1], points[1], points[2]], 0.03
        )
    end
end
