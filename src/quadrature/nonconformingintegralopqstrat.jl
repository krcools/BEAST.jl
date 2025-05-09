struct NonConformingIntegralOpQStrat{S} <: AbstractQuadStrat
    conforming_qstrat::S
end

function quaddata(a, X, Y, tels, bels, qs::NonConformingIntegralOpQStrat)
    return quaddata(a, X, Y, tels, bels, qs.conforming_qstrat)
end

function quadrule(a, 𝒳, 𝒴, i, τ, j, σ, qd,
    qs::NonConformingIntegralOpQStrat)

    if CompScienceMeshes.overlap(τ, σ)
        return NonConformingOverlapQRule(qs.conforming_qstrat)
    end

    for (i,λ) in pairs(faces(τ))
        for (j,μ) in pairs(faces(σ))
            if CompScienceMeshes.overlap(λ, μ)
                return NonConformingTouchQRule(qs.conforming_qstrat, i, j)
    end end end

    # Either positive distance or common vertex, both can
    # be handled directly by the parent quadrature strategy
    return quadrule(a, 𝒳, 𝒴, i, τ, j, σ, qd, qs.conforming_qstrat)    
end