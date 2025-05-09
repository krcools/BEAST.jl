@testitem "norm of constant field" begin
    using CompScienceMeshes
    using LinearAlgebra

    h1 = 0.2
    h2 = 0.176
    m1 = meshcuboid(1.0, 1.0, 1.0, h1; generator=:gmsh)
    m2 = meshcuboid(1.0, 1.0, 1.0, h2; generator=:gmsh)

    E = Maxwell3D.planewave(;direction=point(0,0,1), polarization=point(1,0,0), wavenumber=0.0)
    e = n × E

    cstrat = BEAST.DoubleNumWiltonSauterQStrat{Int64, Int64}(2, 3, 6, 7, 5, 5, 4, 3)
    # cstrat = BEAST.DoubleNumWiltonSauterQStrat{Int64, Int64}(12, 13, 12, 13, 7, 7, 7, 7)
    nstrat = BEAST.NonConformingIntegralOpQStrat(cstrat)

    Id = BEAST.Identity()
    S = Maxwell3D.singlelayer(alpha=1.0, beta=1.0, gamma=1.0)

    X1 = raviartthomas(m1)
    X2 = raviartthomas(m2)

    G11 = assemble(Id, X1, X1)
    G22 = assemble(Id, X2, X2)

    u1 = G11 \ assemble(e, X1)
    u2 = G22 \ assemble(e, X2)

    S11 = assemble(S, X1, X1; quadstrat=cstrat)
    S12 = assemble(S, X1, X2; quadstrat=nstrat)
    S22 = assemble(S, X2, X2; quadstrat=cstrat)

    term1 = real(u1' * S11 * u1)
    term2 = real(u1' * S12 * u2)
    term3 = real(u2' * S22 * u2)

    norm1 = sqrt(term1)
    norm3 =  sqrt(term2)
    testval = norm(term1 - 2*term2 + term3)
    @show testval / norm1
    @test testval / norm1 < 3e-2
end