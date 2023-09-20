# This file contains kernel independent implementations of quadrule for
# the various quadrature strategies defined in quadstrats.jl

function quadrule(op::IntegralOperator, g::RefSpace, f::RefSpace,  i, τ, j, σ, qd,
    qs::DoubleNumSauterQstrat)

    T = eltype(eltype(τ.vertices))
    hits = 0
    dtol = 1.0e3 * eps(T)
    dmin2 = floatmax(T)
    for t in τ.vertices
        for s in σ.vertices
            d2 = LinearAlgebra.norm_sqr(t-s)
            dmin2 = min(dmin2, d2)
            hits += (d2 < dtol)
        end
    end

    hits == 3 && return SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
    hits == 2 && return SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
    hits == 1 && return SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])

    return DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j],)
end


function quadrule(op::IntegralOperator, g::LagrangeRefSpace, f::LagrangeRefSpace,  i, τ, j, σ, qd,
    qs::DoubleNumSauterQstrat)

    T = eltype(eltype(τ.vertices))
    hits = 0
    dtol = 1.0e3 * eps(T)
    dmin2 = floatmax(T)
    for t in τ.vertices
        for s in σ.vertices
            d2 = LinearAlgebra.norm_sqr(t-s)
            dmin2 = min(dmin2, d2)
            hits += (d2 < dtol)
        end
    end

    hits == 3 && return SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
    hits == 2 && return SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
    hits == 1 && return SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])

    return DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j],)
end