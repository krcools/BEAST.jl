@testitem "quadstrat as function" begin
    using CompScienceMeshes
    using LinearAlgebra

    pd = dirname(pathof(BEAST))
    fn = joinpath(pd, "../test/assets/sphere45.in")
    Γ = CompScienceMeshes.readmesh(fn)
    @show length(Γ)

    X = raviartthomas(Γ)
    @show numfunctions(X)
    t = Maxwell3D.singlelayer(wavenumber=1.0)

    qsp(p) = (op, X, Y) -> BEAST.DoubleNumWiltonSauterQStrat(2+p, 3+p, 6+p, 7+p, 5+p, 5+p, 4+p, 3+p)
    Z = map(0:3) do p
        assemble(t, X, X; quadstrat=qsp(p))
    end

    W = map(0:3) do p
        qs = qsp(p)(t, X, X)
        assemble(t, X, X; quadstrat=qs)
    end

    @test all(Z .≈ W)
end

@testitem "quadstrat for linear combinations" begin
    using CompScienceMeshes
    using LinearAlgebra

    pd = dirname(pathof(BEAST))
    fn = joinpath(pd, "../test/assets/sphere45.in")
    Γ = CompScienceMeshes.readmesh(fn)
    @show length(Γ)

    X = raviartthomas(Γ)
    @show numfunctions(X)
    t = Maxwell3D.singlelayer(wavenumber=1.0)
    k = Maxwell3D.doublelayer(wavenumber=1.0)

    qsp(p) = (op, X, Y) -> BEAST.DoubleNumWiltonSauterQStrat(2+p, 3+p, 6+p, 7+p, 5+p, 5+p, 4+p, 3+p)
    
    a = t + k
    qsp0 = qsp(0)
    Z = assemble(a, X, X; quadstrat=qsp0)
    # @which assemble(a, X, X; quadstrat=qsp0)
end