@testitem "operators: interpolators" begin
    using CompScienceMeshes
    using LinearAlgebra

    fn = pkgdir(BEAST, "test/assets", "twoboxes2.geo")
    
    h = 0.25
    m1 = meshgeo(fn, h = h, physical="Gamma1")
    m2 = meshgeo(fn, h = h, physical="Gamma2")
    m3 = meshgeo(fn, h = h, physical="Gamma3")

    CompScienceMeshes.orient(m1)
    CompScienceMeshes.orient(m2)
    CompScienceMeshes.orient(m3)

    m2.vertices = m1.vertices
    m3.vertices = m1.vertices

    Σ = CompScienceMeshes.union(m1, m2, m3)

    ∂Ω0 = CompScienceMeshes.union(-m1, m2)
    ∂Ω1 = CompScienceMeshes.union(m1, -m3)
    ∂Ω2 = CompScienceMeshes.union(-m2, m3)

    U = BEAST.gwpcurl(Σ; order=0, interpolator=true)
    U0 = BEAST.gwpcurl(∂Ω0; order=0, interpolator=true)
    U1 = BEAST.gwpcurl(∂Ω1; order=0, interpolator=true)
    U2 = BEAST.gwpcurl(∂Ω2; order=0, interpolator=true)

    Π = BEAST.InterpolatorOp()
    R0 = assemble(Π, U0, U)
    R1 = assemble(Π, U1, U)
    R2 = assemble(Π, U2, U)

    # Σ is closed, with a junction where the three surfaces meet
    @test length(boundary(Σ)) == 0

    # for the lowest order, interpolation reduces to the signed injection of edges
    for (R, Ui) in ((R0, U0), (R1, U1), (R2, U2))
        @test length(boundary(geometry(Ui))) == 0
        @test size(R) == (numfunctions(Ui), numfunctions(U))
        S = CompScienceMeshes.embedding(Ui.interpolator.edges, U.interpolator.edges)
        @test norm(R - S, Inf) < 1e-10
    end
end