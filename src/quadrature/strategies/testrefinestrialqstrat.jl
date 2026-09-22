struct TestRefinesTrialQStrat{S} <: AbstractQuadStrat
    conforming_qstrat::S
end

function quaddata(a, X, Y, tels, bels, qs::TestRefinesTrialQStrat)
    return quaddata(a, X, Y, tels, bels, qs.conforming_qstrat)
end

function integrate!(a, 𝒳, 𝒴, i, τ, j, σ, qd,
    qs::TestRefinesTrialQStrat,
    out=nothing, test_space=nothing, tptr=nothing, trial_space=nothing, bptr=nothing;
    action::QuadRuleAction)

    hits = _numhits(τ, σ)
    if hits > 0
        qrule = TestRefinesTrialQRule(qs.conforming_qstrat)
        return integrate!(action, out, a, test_space, tptr, τ, trial_space, bptr, σ, qrule)
    end

    return integrate!(a, 𝒳, 𝒴, i, τ, j, σ, qd, qs.conforming_qstrat,
        out, test_space, tptr, trial_space, bptr; action)
end
