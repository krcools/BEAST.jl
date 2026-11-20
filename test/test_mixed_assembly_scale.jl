@testitem "Mixed assembly does not depend on the unit of length" begin
    using CompScienceMeshes
    using CollisionDetection
    using LinearAlgebra
    using SparseArrays
    using StaticArrays

    # A constant kernel makes the separately implemented composed assembly an
    # identity pairing too, with an independent intersection/assembly code path.
    struct MixedAssemblyUnitKernel{T} <: BEAST.Kernel{T} end
    (::MixedAssemblyUnitKernel{T})(x, displacement) where {T} = one(T)
    operators(T) = (Identity(), BEAST.CompSingleKern(dot, MixedAssemblyUnitKernel{T}(), *))
    quadrature(args...) = BEAST.SingleNumQStrat(6)
    gram(op, X, Y) = assemble(op, X, Y; quadstrat=quadrature, verbose=false)

    # Widen comparisons explicitly, while retaining the tolerance of the actual
    # working precision. Never compute a tiny Float16 reference area.
    tolerance(T) = Float64(16eps(T))
    scaled(mesh, s) = Mesh([s * p for p in vertices(mesh)], copy(cells(mesh)))
    function triangle(T, points)
        Mesh(SVector{3,T}.(collect(points)), [CompScienceMeshes.SimplexGraph(1, 2, 3)])
    end

    function check_integral(G, expected, T)
        @test eltype(G) === T
        @test all(isfinite, G)
        @test isfinite(expected)
        if iszero(expected)
            @test iszero(G)
        else
            @test expected > 0
            @test !iszero(G)
            @test isapprox(sum(Float64, G), expected; rtol=tolerance(T), atol=0)
        end
    end

    function nodal_values(space, component, s)
        mesh = geometry(space)
        map(space.fns) do fn
            shape = first(fn)
            node = CompScienceMeshes.indices(mesh, shape.cellid)[shape.refid]
            Float64(vertices(mesh)[node][component]) / Float64(s)
        end
    end

    # Exercise local integration without the octree. Analytic P0 overlap integrals
    # remain meaningful for Float16/Float32 while the tree has a type mismatch.
    function pair_integral(test_mesh, trial_mesh; check_overlap=true)
        T = coordtype(test_mesh)
        X, Y = lagrangecxd0(test_mesh), lagrangecxd0(trial_mesh)
        tcell, bcell = chart(test_mesh, 1), chart(trial_mesh, 1)
        _, tad, _ = BEAST.assemblydata(X)
        _, bad, _ = BEAST.assemblydata(Y)
        trefs, brefs = refspace(X), refspace(Y)
        op = Identity()
        qs = quadrature(op, X, Y)
        qd = BEAST.quaddata(op, trefs, brefs, [tcell], [bcell], qs)
        G = zeros(T, 1, 1)
        if !check_overlap || CompScienceMeshes.overlap(tcell, bcell)
            for cell in CompScienceMeshes.intersection(tcell, bcell)
                BEAST._assemble_restricted_interaction!(
                    op,
                    trefs,
                    brefs,
                    tcell,
                    bcell,
                    cell,
                    tad,
                    bad,
                    1,
                    1,
                    qd,
                    qs,
                    (value, row, col) -> (G[row, col] += value),
                    sqrt(eps(T)),
                )
            end
        end
        G
    end

    @testset "P1 matrix scaling, transpose, and analytic moments" begin
        unit_x = meshrectangle(1.0, 1.0, 0.5)
        unit_y = meshrectangle(1.0, 1.0, 1 / 3)
        X0 = lagrangec0d1(unit_x; dirichlet=false)
        Y0 = lagrangec0d1(unit_y; dirichlet=false)
        reference = Matrix(gram(Identity(), X0, Y0))

        # Preserve connectivity and function order when comparing matrix entries.
        for s in (1.0, 1e-3, 1e-6, 1e3), op in operators(Float64)
            @testset "$(typeof(op)), side=$s" begin
                X = lagrangec0d1(scaled(unit_x, s); dirichlet=false)
                Y = lagrangec0d1(scaled(unit_y, s); dirichlet=false)
                @test geometry(X) != geometry(Y)
                G = gram(op, X, Y)
                check_integral(G, s^2, Float64)
                normalized = Matrix(G) / s^2
                @test normalized ≈ reference rtol=tolerance(Float64) atol=0
                @test G ≈ transpose(gram(op, Y, X)) rtol=tolerance(Float64) atol=0

                # Independent moments detect wrong entries even when sums agree.
                u = nodal_values(X, 1, s)
                v = nodal_values(Y, 1, s)
                w = nodal_values(Y, 2, s)
                @test sum(transpose(u) * normalized) ≈ 1 / 2 rtol=tolerance(Float64)
                @test transpose(u) * normalized * v ≈ 1 / 3 rtol=tolerance(Float64)
                @test transpose(u) * normalized * w ≈ 1 / 4 rtol=tolerance(Float64)
            end
        end
    end

    @testset "Vector basis scaling" begin
        unit_x = meshrectangle(1.0, 1.0, 0.5)
        unit_y = meshrectangle(1.0, 1.0, 1 / 3)
        for basis in (BEAST.nedelec, raviartthomas)
            reference = gram(Identity(), basis(unit_x), basis(unit_y))
            @test !iszero(reference)
            for s in (1e-3, 1e-6), op in operators(Float64)
                X, Y = basis(scaled(unit_x, s)), basis(scaled(unit_y, s))
                G = gram(op, X, Y)
                @test all(isfinite, nonzeros(G))
                # Unit edge/flux bases scale as 1/s, cancelling the area factor.
                @test G ≈ reference rtol=tolerance(Float64) atol=0
                @test G ≈ transpose(gram(op, Y, X)) rtol=tolerance(Float64) atol=0
            end
        end
    end

    # Reordered coincident triangles force the mixed dispatcher: the point sets
    # coincide, but the indexed meshes do not compare equal.
    geometry_cases = (
        (:coincident, ((1, 0, 0), (0, 1, 0), (0, 0, 0)), 1 / 2),
        (:reversed, ((0, 1, 0), (1, 0, 0), (0, 0, 0)), 1 / 2),
        (:contained, ((1 / 4, 1 / 4, 0), (1 / 2, 1 / 4, 0), (1 / 4, 1 / 2, 0)), 1 / 32),
        (:partial, ((1 / 2, 0, 0), (3 / 2, 0, 0), (1 / 2, 1, 0)), 1 / 8),
        (:edge_touch, ((0, 0, 0), (1, 0, 0), (0, -1, 0)), 0),
        (:vertex_touch, ((0, 0, 0), (-1, 0, 0), (0, -1, 0)), 0),
        (:disjoint, ((2, 0, 0), (3, 0, 0), (2, 1, 0)), 0),
        (:different_planes, ((0, 0, 1), (1, 0, 1), (0, 1, 1)), 0),
    )

    @testset "Typed local integration and geometric corner cases" begin
        for T in (Float16, Float32, Float64), s in (T(1 / 2), one(T), T(2))
            @testset "$T, side=$s" begin
                a = scaled(triangle(T, ((0, 0, 0), (1, 0, 0), (0, 1, 0))), s)
                @test coordtype(a) === T
                for (name, points, area) in geometry_cases
                    @testset "$name" begin
                        b = scaled(triangle(T, points), s)
                        @test coordtype(b) === T
                        expected = Float64(s)^2 * area
                        G = pair_integral(a, b)
                        check_integral(G, expected, T)
                        @test G ≈ transpose(pair_integral(b, a)) rtol=tolerance(T) atol=0
                    end
                end
            end
        end
    end

    @testset "Relative cutoff uses the smaller parent" begin
        for T in (Float16, Float32, Float64)
            # h^2 < sqrt(eps(T)), but the entire contained cell must be kept.
            # Using max instead of min would discard it. Coordinates are dyadic.
            h = ldexp(one(T), -cld(precision(T), 4))
            a = triangle(T, ((0, 0, 0), (1, 0, 0), (0, 1, 0)))
            b = triangle(
                T,
                (
                    (T(1 / 4), T(1 / 4), 0),
                    (T(1 / 4) + h, T(1 / 4), 0),
                    (T(1 / 4), T(1 / 4) + h, 0),
                ),
            )
            expected = Float64(h)^2 / 2
            @test expected > 0
            @test h^2 < sqrt(eps(T))
            # Containment is known analytically. Isolate the assembly cutoff
            # from the upstream overlap predicate's separate Float16 limitation.
            check_integral(pair_integral(a, b; check_overlap=false), expected, T)
            check_integral(pair_integral(b, a; check_overlap=false), expected, T)
            @test CompScienceMeshes.overlap(chart(a, 1), chart(b, 1))
            if T === Float16
                # The separating-axis tolerance is currently set by the larger
                # triangle and wrongly rejects this pair in the opposite order.
                @test_broken CompScienceMeshes.overlap(chart(b, 1), chart(a, 1))
            else
                @test CompScienceMeshes.overlap(chart(b, 1), chart(a, 1))
            end

            # On equal parents, exercise the intended sliver-filtering policy
            # below the cutoff and an intersection comfortably above it.
            for factor in (1, 4)
                delta = T(factor) * h
                offset = one(T) - delta
                partial = triangle(
                    T, ((offset, 0, 0), (offset + one(T), 0, 0), (offset, 1, 0))
                )
                expected = factor == 1 ? 0.0 : Float64(delta)^2 / 2
                check_integral(pair_integral(a, partial), expected, T)
                check_integral(pair_integral(partial, a), expected, T)
            end
        end
    end

    @testset "End-to-end P0 geometry cases in both assembly paths" begin
        for s in (1.0, 1e-3, 1e-6), op in operators(Float64)
            a = scaled(triangle(Float64, ((0, 0, 0), (1, 0, 0), (0, 1, 0))), s)
            for (name, points, area) in geometry_cases
                @testset "$(typeof(op)), $name, side=$s" begin
                    b = scaled(triangle(Float64, points), s)
                    X, Y = lagrangecxd0(a), lagrangecxd0(b)
                    @test geometry(X) != geometry(Y)
                    G = gram(op, X, Y)
                    check_integral(G, s^2 * area, Float64)
                    @test G ≈ transpose(gram(op, Y, X)) rtol=tolerance(Float64) atol=0
                end
            end
        end
    end

    @testset "Known lower-precision octree constructor limitation" begin
        for T in (Float16, Float32), op in operators(T)
            @testset "$(typeof(op)), $T" begin
                X = lagrangec0d1(meshrectangle(one(T), one(T), T(1 / 2)); dirichlet=false)
                Y = lagrangec0d1(
                    meshrectangle(one(T), one(T), one(T) / T(3)); dirichlet=false
                )
                @test coordtype(geometry(X)) === T
                @test scalartype(X) === T
                G = try
                    gram(op, X, Y)
                catch err
                    # The resolved CollisionDetection constructor requires a
                    # Float64 expansion_ratio. Unrelated errors must still fail.
                    err isa MethodError && err.f === CollisionDetection.Octree ||
                        rethrow()
                    nothing
                end
                # A successful assembly must cause an unexpected pass, so this
                # marker is removed when the constructor mismatch is fixed.
                @test_broken G !== nothing
                if G !== nothing
                    check_integral(G, 1.0, T)
                end
            end
        end
    end

    @testset "Known absolute cutoff in triangle intersection" begin
        s = 1e-8
        a = scaled(triangle(Float64, ((0, 0, 0), (1, 0, 0), (0, 1, 0))), s)
        b = scaled(triangle(Float64, ((1 / 2, 0, 0), (3 / 2, 0, 0), (1 / 2, 1, 0))), s)
        expected = s^2 / 8
        @test isfinite(expected) && expected > 0
        for op in operators(Float64)
            G = gram(op, lagrangecxd0(a), lagrangecxd0(b))
            @test all(isfinite, nonzeros(G))
            # intersection/intersection2 still filter using an absolute eps().
            # A dependency fix must turn this marker into an unexpected pass.
            @test_broken isapprox(
                sum(Float64, G), expected; rtol=tolerance(Float64), atol=0
            )
        end
    end
end
