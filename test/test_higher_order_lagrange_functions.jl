@testitem "lagrangecx order=3 - global" begin
    using CompScienceMeshes

    projectdir = joinpath(dirname(pathof(BEAST)),"..")
    m = readmesh(joinpath(projectdir, "test/assets/sphere2.in"))

    lagspace3 = BEAST.lagrangecx(m, order=3)
    @test numfunctions(lagspace3) == 10 * length(m)
    @test refspace(lagspace3) == BEAST.LagrangeRefSpace{Float64,3,3,10}()
end

@testitem "lagrangecxd3: local, ref" begin
    using CompScienceMeshes

    T = Float64
    s = simplex(
        point(1,0,0),
        point(0,1,0),
        point(0,0,0),
    )
    p = center(s)

    ϕ = BEAST.LagrangeRefSpace{Float64,3,3,10}()
    A = ϕ(p)

    # test the dimension
    @test length(A) == binomial(5,2)

    # test the partition of unity property
    valp = sum(a.value for a in A)
    crlp = sum(a.curl for a in A)
    @test valp ≈ 1
    @test crlp ≈ point(0,0,0) atol=sqrt(eps(T))
    @show crlp
end


@testitem "lagrangecxd3: local, generic simplex" begin
    using CompScienceMeshes

    s = simplex(
        point(3,0,0),
        point(0,2,1),
        point(-1,-1,-1),
    )
    p = center(s)

    ϕ = BEAST.LagrangeRefSpace{Float64,3,3,10}()
    A = ϕ(p)

    # test the dimension
    @test length(A) == binomial(5,2)

    # test the partition of unity property
    valp = sum(a.value for a in A)
    @test valp ≈ 1
end

@testitem "lagrangecxd3: self-interpolate" begin
    using CompScienceMeshes
    using LinearAlgebra

    ϕ = BEAST.LagrangeRefSpace{Float64,3,3,10}()
    ch = simplex(
        point(3,0,0),
        point(0,2,1),
        point(-1,-1,-1))
    ch1toch2 = simplex(
        point(1,0),
        point(0,1),
        point(0,0))

    Q = BEAST.interpolate(ϕ, ch, ϕ, ch, ch1toch2)
    # display(round.(Q, digits=3))
    Id = Matrix{Float64}(I, 10, 10)
    @test Q ≈ Id
end


@testitem "lagrnacexcd3: interpolate generic poly" begin

    using CompScienceMeshes
    using LinearAlgebra

    ϕ = BEAST.LagrangeRefSpace{Float64,3,3,10}()
    ch = simplex(
        point(3,0,0),
        point(0,2,1),
        point(-1,-1,-1))

    function fields(p)
        u, v = parametric(p)
        return [
            1 + u + v + u^2 + u*v + v^2 + u^3 + u^2*v + u*v^2 + v^3,
            1 + u + u^2 + u^3,
            1 + u + u^2,
            1 + u,
            1
        ]
    end

    Q = BEAST.interpolate(fields, ϕ, ch)
    nbds = [
        neighborhood(ch, (0.2, 0.2)),
        neighborhood(ch, (0.2, 0.6)),
        neighborhood(ch, (0.6, 0.2)),
        neighborhood(ch, (1/3, 1/3)),
        neighborhood(ch, (0.0, 0.0))]

    for p in nbds
        basis = ϕ(p)

        vals = fields(p)
        for j in eachindex(vals)
            val1 = vals[j]
            val2 = sum(Q[j,i] * b.value for (i,b) in enumerate(basis))
            @show val1 val2
        end
        println()
    end
end