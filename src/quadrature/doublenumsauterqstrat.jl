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
function quaddata(op::IntegralOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts::Vector{<:CompScienceMeshes.Simplex{<:Any, 3}}, trial_charts::Vector{<:CompScienceMeshes.Simplex{<:Any, 3}}, qs::DoubleNumSauterQstrat)

    #The combinations of rules (6,7) and (5,7 are) BAAAADDDD
    # they result in many near singularity evaluations with any
    # resemblence of accuracy going down the drain! Simply don't!
    # (same for (5,7) btw...).
    t_qp = quadpoints(test_local_space,  test_charts,  (qs.outer_rule,))
    b_qp = quadpoints(trial_local_space, trial_charts, (qs.inner_rule,))

    sing_qp = (SauterSchwab3D._legendre(qs.sauter_schwab_common_vert,0,1), 
               SauterSchwab3D._shunnham2D(qs.sauter_schwab_common_edge),
               SauterSchwab3D._shunnham3D(qs.sauter_schwab_common_face),
               SauterSchwab3D._shunnham4D(qs.sauter_schwab_common_tetr),)
    return (tpoints=t_qp, bpoints=b_qp, sing_qp=sing_qp)
end
function quaddata(op::IntegralOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts::Vector{<:CompScienceMeshes.Simplex{<:Any, 2}}, trial_charts::Vector{<:CompScienceMeshes.Simplex{<:Any, 3}}, qs::DoubleNumSauterQstrat)

    #The combinations of rules (6,7) and (5,7 are) BAAAADDDD
    # they result in many near singularity evaluations with any
    # resemblence of accuracy going down the drain! Simply don't!
    # (same for (5,7) btw...).
    t_qp = quadpoints(test_local_space,  test_charts,  (qs.outer_rule,))
    b_qp = quadpoints(trial_local_space, trial_charts, (qs.inner_rule,))

    sing_qp = (SauterSchwab3D._legendre(qs.sauter_schwab_common_vert,0,1), 
               SauterSchwab3D._shunnham2D(qs.sauter_schwab_common_edge),
               SauterSchwab3D._shunnham3D(qs.sauter_schwab_common_face),
               SauterSchwab3D._shunnham4D(qs.sauter_schwab_common_tetr),)
    return (tpoints=t_qp, bpoints=b_qp, sing_qp=sing_qp)
end

function quaddata(op::IntegralOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts::Vector{<:CompScienceMeshes.Simplex{<:Any, 3}}, trial_charts::Vector{<:CompScienceMeshes.Simplex{<:Any, 2}}, qs::DoubleNumSauterQstrat)

    #The combinations of rules (6,7) and (5,7 are) BAAAADDDD
    # they result in many near singularity evaluations with any
    # resemblence of accuracy going down the drain! Simply don't!
    # (same for (5,7) btw...).
    t_qp = quadpoints(test_local_space,  test_charts,  (qs.outer_rule,))
    b_qp = quadpoints(trial_local_space, trial_charts, (qs.inner_rule,))

    sing_qp = (SauterSchwab3D._legendre(qs.sauter_schwab_common_vert,0,1), 
               SauterSchwab3D._shunnham2D(qs.sauter_schwab_common_edge),
               SauterSchwab3D._shunnham3D(qs.sauter_schwab_common_face),
               SauterSchwab3D._shunnham4D(qs.sauter_schwab_common_tetr),)
    return (tpoints=t_qp, bpoints=b_qp, sing_qp=sing_qp)
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

struct _TransposedStrat{A}
    strat::A
end 

function quadrule(op::IntegralOperator, g::RefSpace, f::RefSpace, 
    i, τ::CompScienceMeshes.Simplex{<:Any, 3},
    j, σ::CompScienceMeshes.Simplex{<:Any, 3}, 
    qd, qs::DoubleNumSauterQstrat) 
    qr_volume(op, g, f, i, τ, j, σ, qd, qs)
end
function quadrule(op::IntegralOperator, g::RefSpace, f::RefSpace, 
    i, τ::CompScienceMeshes.Simplex{<:Any, 3},
    j, σ::CompScienceMeshes.Simplex{<:Any, 2}, 
    qd, qs::DoubleNumSauterQstrat) 
    qr_boundary(op, g, f, i, τ, j, σ, qd, qs)
end

function quadrule(op::IntegralOperator, g::RefSpace, f::RefSpace, 
    i, τ::CompScienceMeshes.Simplex{<:Any, 2},
    j, σ::CompScienceMeshes.Simplex{<:Any, 3}, 
    qd, qs::DoubleNumSauterQstrat) 
    _TransposedStrat(qr_boundary(op, g, f, i, τ, j, σ, qd, qs))
end

function qr_volume(op::IntegralOperator, g::RefSpace, f::RefSpace, i, τ, j, σ, qd, qs)

    dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))

    hits = 0
    idx_t = Int64[]
    idx_s = Int64[]
    sizehint!(idx_t,4)
    sizehint!(idx_s,4)
    dmin2 = floatmax(eltype(eltype(τ.vertices)))
    D = dimension(τ)+dimension(σ)
    for (i,t) in enumerate(τ.vertices)
        for (j,s) in enumerate(σ.vertices)
            d2 = LinearAlgebra.norm_sqr(t-s)
            d = norm(t-s)
            dmin2 = min(dmin2, d2)
            # if d2 < dtol
            if d < dtol
                push!(idx_t,i)
                push!(idx_s,j)
                hits +=1
                break
            end
        end
    end

    #singData = SauterSchwab3D.Singularity{D,hits}(idx_t, idx_s )
   @assert hits <= 4

    hits == 4 && return SauterSchwab3D.CommonVolume6D_S(SauterSchwab3D.Singularity6DVolume(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[4]))
    hits == 3 && return SauterSchwab3D.CommonFace6D_S(SauterSchwab3D.Singularity6DFace(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
    hits == 2 && return SauterSchwab3D.CommonEdge6D_S(SauterSchwab3D.Singularity6DEdge(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3],qd.sing_qp[4]))
    hits == 1 && return SauterSchwab3D.CommonVertex6D_S(SauterSchwab3D.Singularity6DPoint(idx_t,idx_s),qd.sing_qp[3])



    return DoubleQuadRule(
        qd[1][1,i],
        qd[2][1,j])

end


function qr_boundary(op::IntegralOperator, g::RefSpace, f::RefSpace, i, τ, j,  σ, qd, qs)

    dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))

    hits = 0
    idx_t = Int64[]
    idx_s = Int64[]
    sizehint!(idx_t,4)
    sizehint!(idx_s,4)
    dmin2 = floatmax(eltype(eltype(τ.vertices)))
    D = dimension(τ)+dimension(σ)
    for (i,t) in enumerate(τ.vertices)
        for (j,s) in enumerate(σ.vertices)
            d2 = LinearAlgebra.norm_sqr(t-s)
            d = norm(t-s)
            dmin2 = min(dmin2, d2)
            # if d2 < dtol
            if d < dtol
                push!(idx_t,i)
                push!(idx_s,j)
                hits +=1
                break
            end
        end
    end

    @assert hits <= 3
    #singData = SauterSchwab3D.Singularity{D,hits}(idx_t, idx_s )
   

    hits == 3 && return SauterSchwab3D.CommonFace5D_S(SauterSchwab3D.Singularity5DFace(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
    hits == 2 && return SauterSchwab3D.CommonEdge5D_S(SauterSchwab3D.Singularity5DEdge(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
    hits == 1 && return SauterSchwab3D.CommonVertex5D_S(SauterSchwab3D.Singularity5DPoint(idx_t,idx_s),(qd.sing_qp[3],qd.sing_qp[2]))


    return DoubleQuadRule(
        qd[1][1,i],
        qd[2][1,j])

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