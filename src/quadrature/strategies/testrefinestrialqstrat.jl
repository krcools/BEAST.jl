struct TestRefinesTrialQStrat{S} <: AbstractQuadStrat
    conforming_qstrat::S
end

function quaddata(a, X, Y, tels, bels, qs::TestRefinesTrialQStrat)
    return quaddata(a, X, Y, tels, bels, qs.conforming_qstrat)
end

function quadrule(a, 𝒳, 𝒴, i, τ, j, σ, qd,
    qs::TestRefinesTrialQStrat)

    hits = _numhits(τ, σ)
    if hits > 0
        return TestRefinesTrialQRule(qs.conforming_qstrat)
    end

    return quadrule(a, 𝒳, 𝒴, i, τ, j, σ, qd, qs.conforming_qstrat)    
end