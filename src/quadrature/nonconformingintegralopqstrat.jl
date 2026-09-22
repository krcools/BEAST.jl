struct NonConformingIntegralOpQStrat{S} <: AbstractQuadStrat
    conforming_qstrat::S
end

function quaddata(a, X, Y, tels, bels, qs::NonConformingIntegralOpQStrat)
    return quaddata(a, X, Y, tels, bels, qs.conforming_qstrat)
end

function integrate!(a, 𝒳, 𝒴, i, τ, j, σ, qd,
    qs::NonConformingIntegralOpQStrat,
    out=nothing, test_space=nothing, tptr=nothing, trial_space=nothing, bptr=nothing;
    action::QuadRuleAction)

    if CompScienceMeshes.overlap(τ, σ)
        qrule = NonConformingOverlapQRule(qs.conforming_qstrat)
        return integrate!(action, out, a, test_space, tptr, τ, trial_space, bptr, σ, qrule)
    end

    for (i,λ) in pairs(faces(τ))
        for (j,μ) in pairs(faces(σ))
            if CompScienceMeshes.overlap(λ, μ)
                qrule = NonConformingTouchQRule(qs.conforming_qstrat, i, j)
                return integrate!(action, out, a, test_space, tptr, τ, trial_space, bptr, σ, qrule)
    end end end

    # Either positive distance or common vertex, both can
    # be handled directly by the parent quadrature strategy
    return integrate!(a, 𝒳, 𝒴, i, τ, j, σ, qd, qs.conforming_qstrat,
        out, test_space, tptr, trial_space, bptr; action)
end
