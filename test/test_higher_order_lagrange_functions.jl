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