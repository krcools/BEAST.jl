struct TestRefinesTrialQStrat{S}
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

    # if CompScienceMeshes.overlap(τ, σ)
    #     return NonConformingOverlapQRule(qs.conforming_qstrat)
    # end

    # for (i,λ) in pairs(faces(τ))
    #     for (j,μ) in pairs(faces(σ))
    #         if CompScienceMeshes.overlap(λ, μ)
    #             return NonConformingOverlapQRule(qs.conforming_qstrat)
    # end end end

    return quadrule(a, 𝒳, 𝒴, i, τ, j, σ, qd, qs.conforming_qstrat)    
end