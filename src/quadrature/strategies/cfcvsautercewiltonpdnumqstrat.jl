struct CommonFaceVertexSauterCommonEdgeWiltonPostitiveDistanceNumQStrat{R,S} <: AbstractQuadStrat
    outer_rule_far::R
    inner_rule_far::R
    outer_rule_near::R
    inner_rule_near::R
    sauter_schwab_common_tetr::S
    sauter_schwab_common_face::S
    sauter_schwab_common_edge::S
    sauter_schwab_common_vert::S
end

function quaddata(op::IntegralOperator, test_local_space, bsis_local_space,
    test_charts, bsis_charts, qs::CommonFaceVertexSauterCommonEdgeWiltonPostitiveDistanceNumQStrat)

    return quaddata(op, test_local_space, bsis_local_space,
        test_charts, bsis_charts, DoubleNumWiltonSauterQStrat(
            qs.outer_rule_far,
            qs.inner_rule_far,
            qs.outer_rule_near,
            qs.inner_rule_near,
            qs.sauter_schwab_common_tetr,
            qs.sauter_schwab_common_face,
            qs.sauter_schwab_common_edge,
            qs.sauter_schwab_common_vert,))
end

function integrate!(op::IntegralOperator, g, f,  i, τ, j, σ,
    qd, qs::CommonFaceVertexSauterCommonEdgeWiltonPostitiveDistanceNumQStrat,
    out=nothing, test_space=nothing, tptr=nothing, trial_space=nothing, bptr=nothing;
    action::QuadRuleAction)

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

    if hits == 3
        qrule = SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
        return integrate!(action, out, op, test_space, tptr, τ, trial_space, bptr, σ, qrule)
    end
    if hits == 2
        qrule = WiltonSERule(
            qd.tpoints[2,i],
            SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2]),)
        return integrate!(action, out, op, test_space, tptr, τ, trial_space, bptr, σ, qrule)
    end
    if hits == 1
        qrule = SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])
        return integrate!(action, out, op, test_space, tptr, τ, trial_space, bptr, σ, qrule)
    end

    h2 = volume(σ)
    xtol2 = 0.2 * 0.2
    k2 = abs2(gamma(op))
    if max(dmin2*k2, dmin2/16h2) < xtol2
        qrule = WiltonSERule(
            qd.tpoints[2,i],
            DoubleQuadRule(
                qd.tpoints[2,i],
                qd.bpoints[2,j],),)
        return integrate!(action, out, op, test_space, tptr, τ, trial_space, bptr, σ, qrule)
    end

    qrule = DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j],)
    return integrate!(action, out, op, test_space, tptr, τ, trial_space, bptr, σ, qrule)
end
