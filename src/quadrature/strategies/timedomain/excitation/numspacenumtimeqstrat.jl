struct NumSpaceNumTimeQStrat{S,T} <: AbstractQuadStrat
    space_rule::S
    time_rule::T
end


function quaddata(exc::TDFunctional,
    testrefs, timerefs,
    testels, timeels, quadstrat::NumSpaceNumTimeQStrat)

    r = quadstrat.space_rule
    s = quadstrat.time_rule

    testqd = quadpoints(testrefs, testels, (r,))
    timeqd = quadpoints(timerefs, timeels, (s,))

    testqd, timeqd

end


function quadrule(exc::TDFunctional,
    testrefs, timerefs,
    p, τ, r, ρ,
    qd, quadstrat::NumSpaceNumTimeQStrat)

    MultiQuadStrategy(
        qd[1][1,p],
        SingleQuadStrategy2(
            qd[2][1,r]
        )
    )

end


function quaddata(excitation::TDFunctional,
    test_refspace, time_refspace::DiracBoundary,
    test_elements, time_elements, quadstrat::NumSpaceNumTimeQStrat)

    r = quadstrat.space_rule
    test_quad_data = quadpoints(test_refspace, test_elements, (r,))

    test_quad_data, nothing
end


function quadrule(exc::TDFunctional,
    testrefs, timerefs::DiracBoundary,
    p, τ, r, ρ,
    qd, quadstrat::NumSpaceNumTimeQStrat)

    MultiQuadStrategy(
        qd[1][1,p],
        nothing
    )

end


@testitem "tdexc: multiquadqrule" begin
    using CompScienceMeshes
    using LinearAlgebra

    fn = joinpath(dirname(pathof(BEAST)), "../test/assets/sphere45.in")

    mesh = readmesh(fn)
    RT = raviartthomas(mesh)

    Δt = 0.1
    Nt = 200
    T = timebasisshiftedlagrange(Δt, Nt, 3)
    U = timebasisdelta(Δt, Nt)

    X = RT ⊗ U

    duration = 2 * 20 * Δt
    delay = 1.5 * duration
    amplitude = 1.0
    gaussian = creategaussian(duration, delay, amplitude)
    direction, polarisation = ẑ, x̂
    E = planewave(polarisation, direction, derive(gaussian), 1.0)

    qs1 = BEAST.NumSpaceNumTimeQStrat(2, 10)
    qs2 = BEAST.NumSpaceNumTimeQStrat(6, 20)

    b1 = assemble(E, X; quadstrat=qs1)
    b2 = assemble(E, X; quadstrat=qs2)

    @test norm(b1-b2, Inf) < 1e-4
end