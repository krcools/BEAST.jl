struct SelfSauterOtherwiseDNumQStrat{R,S} <: AbstractQuadStrat
    outer_rule::R
    inner_rule::R
    sauter_schwab_common::S
end

function quaddata(op::IntegralOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts, trial_charts, qs::SelfSauterOtherwiseDNumQStrat)

    T = coordtype(test_charts[1])

    tqd = quadpoints(test_local_space,  test_charts,  (qs.outer_rule,))
    bqd = quadpoints(trial_local_space, trial_charts, (qs.inner_rule,))
     
    leg = (
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common,0,1)),
      )

    return (tpoints=tqd, bpoints=bqd, gausslegendre=leg)
end


function quadrule(op::IntegralOperator, g::RefSpace, f::RefSpace,  i, τ, j, σ, qd,
    qs::SelfSauterOtherwiseDNumQStrat)

    T = eltype(eltype(τ.vertices))
    hits = 0
    dtol = 1.0e3 * eps(T)
    dmin2 = floatmax(T)
    for t in τ.vertices
        for s in σ.vertices
            d2 = LinearAlgebra.norm_sqr(t-s)
            d = norm(t-s)
            dmin2 = min(dmin2, d2)
            # hits += (d2 < dtol)
            hits += (d < dtol)
        end
    end

    @assert hits <= 3

    hits == 3 && return SauterSchwabQuadrature.CommonFace(qd.gausslegendre[1])

    return DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j],)
end