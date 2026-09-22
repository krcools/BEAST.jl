struct NonConfTestBaryRefOfTrialQStrat{P} <: BEAST.AbstractQuadStrat
    conforming_qstrat::P
end

function BEAST.quaddata(a, X, Y, tels, bels, qs::NonConfTestBaryRefOfTrialQStrat)
    return BEAST.quaddata(a, X, Y, tels, bels, qs.conforming_qstrat)
end

function BEAST.integrate!(a, 𝒳, 𝒴, i, τ, j, σ, qd,
    quadstrat::NonConfTestBaryRefOfTrialQStrat,
    out=nothing, test_space=nothing, tptr=nothing, trial_space=nothing, bptr=nothing;
    action::BEAST.QuadRuleAction)

    # return TestInBaryRefOfTrialQRule(quadstrat.conforming_qstrat)
    nh = BEAST._numhits(τ, σ)
    if nh > 0
        qrule = TestInBaryRefOfTrialQRule(quadstrat.conforming_qstrat)
        return BEAST.integrate!(action, out, a, test_space, tptr, τ, trial_space, bptr, σ, qrule)
    end
    return BEAST.integrate!(a, 𝒳, 𝒴, i, τ, j, σ, qd,
        quadstrat.conforming_qstrat,
        out, test_space, tptr, trial_space, bptr; action)
end

@testitem "NonConfTestBaryRefOfTrialQStrat" begin
    using BEAST, Test
    using CompScienceMeshes, LinearAlgebra

    # @show pathof(CompScienceMeshes)

    fnm = joinpath(dirname(pathof(BEAST)), "../test/assets/sphere45.in")
    Γ1 = BEAST.readmesh(fnm)
    Γ2 = deepcopy(Γ1)

    X = raviartthomas(Γ1)
    Y1 = buffachristiansen(Γ1)
    Y = buffachristiansen(Γ2)

    geoX = BEAST.geometry(X)
    geoY = BEAST.geometry(Y)
    geoY1 = BEAST.geometry(Y1)

    @test geoX === geoY1.parent
    @test CompScienceMeshes.refines(geoY1, geoX)
    @test !CompScienceMeshes.refines(geoY, geoX)
    # @test geoX === CompScienceMeshes.parent(geoY1)
    # @test geoX !== CompScienceMeshes.parent(geoY)

    K = Maxwell3D.doublelayer(gamma=1.0)
    qs1 = BEAST.DoubleNumWiltonSauterQStrat(2, 3, 6, 7, 5, 5, 4, 3)
    qs2 = BEAST.NonConformingIntegralOpQStrat(qs1)
    qs3 = BEAST.NonConfTestBaryRefOfTrialQStrat(qs1)

    @time Kyx1 = assemble(K, Y, X; quadstrat=qs1)
    @time Kyx2 = assemble(K, Y, X; quadstrat=qs2)
    @time Kyx3 = assemble(K, Y, X; quadstrat=qs3)
    @time Kyx4 = assemble(K, Y1, X; quadstrat=qs1)
    
    @test norm(Kyx3 - Kyx4) < 1e-12
    @test norm(Kyx1 - Kyx2) < 0.05
    @test norm(Kyx1 - Kyx3) < 0.05
    @test norm(Kyx2 - Kyx3) < 0.002
    @test norm(Kyx2 - Kyx4) < 0.002
end
