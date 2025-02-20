struct DoubleNumSauterQstrat{R,S} <: AbstractQuadStrat
    outer_rule::R
    inner_rule::R
    sauter_schwab_common_tetr::S
    sauter_schwab_common_face::S
    sauter_schwab_common_edge::S
    sauter_schwab_common_vert::S
end

function quaddata(op::IntegralOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts, trial_charts, qs::DoubleNumSauterQstrat)

    T = coordtype(test_charts[1])

    tqd = quadpoints(test_local_space,  test_charts,  (qs.outer_rule,))
    bqd = quadpoints(trial_local_space, trial_charts, (qs.inner_rule,))
     
    leg = (
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_vert,0,1)),
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_edge,0,1)),
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_face,0,1)),
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_tetr,0,1)),
      )

    return (tpoints=tqd, bpoints=bqd, gausslegendre=leg)
end


function quadrule(op::IntegralOperator, g::RefSpace, f::RefSpace,
    i, τ::CompScienceMeshes.Simplex{<:Any, 2},
    j, σ::CompScienceMeshes.Simplex{<:Any, 2},
    qd, qs::DoubleNumSauterQstrat)

    hits = _numhits(τ, σ)
    @assert hits <= 3

    hits == 3 && return SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
    hits == 2 && return SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
    hits == 1 && return SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])

    return DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j],)
end

function quadrule(op::IntegralOperator, g::RefSpace, f::RefSpace,
    i, τ::CompScienceMeshes.Quadrilateral,
    j, σ::CompScienceMeshes.Quadrilateral,
    qd, qs::DoubleNumSauterQstrat)

    hits = _numhits(τ, σ)
    @assert hits != 3
    @assert hits <= 4

    hits == 4 && return SauterSchwabQuadrature.CommonFaceQuad(qd.gausslegendre[3])
    hits == 2 && return SauterSchwabQuadrature.CommonEdgeQuad(qd.gausslegendre[2])
    hits == 1 && return SauterSchwabQuadrature.CommonVertexQuad(qd.gausslegendre[1])

    return DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j],)
end


function _numhits(τ, σ)
    T = coordtype(τ)
    hits = 0
    dtol = 1.0e3 * eps(T)
    dmin2 = floatmax(T)
    for t in vertices(τ)
        for s in vertices(σ)
            d2 = LinearAlgebra.norm_sqr(t-s)
            d = norm(t-s)
            dmin2 = min(dmin2, d2)
            hits += (d < dtol)
        end
    end
    return hits
end