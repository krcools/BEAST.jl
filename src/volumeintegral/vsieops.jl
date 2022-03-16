abstract type VSIEOperator <: IntegralOperator end
abstract type VolumeSurfaceOperator <: VSIEOperator end
abstract type BoundarySurfaceOperator <: VSIEOperator end
abstract type VSIEOperatorT <: VSIEOperator end
abstract type VolumeSurfaceOperatorT <: VolumeSurfaceOperator end
abstract type BoundarySurfaceOperatorT <: BoundarySurfaceOperator end

struct KernelValsVSIE{T,U,P,Q,K}
    gamma::U
    vect::P
    dist::T
    green::U
    gradgreen::Q
    tau::K
end

function kernelvals(viop::VSIEOperator, p ,q)
    Y = viop.gamma;
    r = cartesian(p)-cartesian(q)
    R = norm(r)
    yR = Y*R

    expn = exp(-yR)
    green = expn / (4*pi*R)
    gradgreen = - (Y +1/R)*green/R*r

    tau = viop.tau(cartesian(q))

    KernelValsVSIE(Y,r,R, green, gradgreen,tau)
end

struct VSIESingleLayer{T,U,P} <: VolumeSurfaceOperator
    gamma::T
    α::U
    β::U
    tau::P
end

struct VSIEBoundary{T,U,P} <: BoundarySurfaceOperator
    gamma::T
    α::U
    tau::P
end

struct VSIEDoubleLayer{T,U,P} <: VolumeSurfaceOperator
    gamma::T
    α::U
    tau::P
end

#=
struct VSIESingleLayer2{T,U,P} <: VolumeSurfaceOperator
    gamma::T
    α::U
    β::U
    tau::P
end

struct VSIEBoundaryT{T,U,P} <: BoundarySurfaceOperatorT
    gamma::T
    α::U
    tau::P
end

struct VSIEDoubleLayerT{T,U,P} <: VolumeSurfaceOperatorT
    gamma::T
    α::U
    tau::P
end
=#
scalartype(op::VSIEOperator) = typeof(op.gamma)

export VSIE

struct VSIEIntegrand{S,T,O,K,L}
    test_tetrahedron_element::S
    trial_triangle_element::T
    op::O
    test_local_space::K
    trial_local_space::L
end

#=
struct VSIEIntegrandT{S,T,O,K,L}
    test_tetrahedron_element::S
    trial_triangle_element::T
    op::O
    test_local_space::K
    trial_local_space::L
end
=#

function (igd::VSIEIntegrand)(u,v)

    #mesh points
    tgeo = neighborhood(igd.test_tetrahedron_element,v)
    bgeo = neighborhood(igd.trial_triangle_element,u)

    #kernel values
    kerneldata = kernelvals(igd.op,tgeo,bgeo)

    #values & grad/div/curl of local shape functions
    tval = igd.test_local_space(tgeo) 
    bval = igd.trial_local_space(bgeo)

    #jacobian
    j = jacobian(tgeo) * jacobian(bgeo)
    
    integrand(igd.op, kerneldata,tval,tgeo,bval,tgeo) * j
end


#=
function (igd::VSIEIntegrandT)(u,v)

    #mesh points
    tgeo = neighborhood(igd.test_tetrahedron_element,u)
    bgeo = neighborhood(igd.trial_triangle_element,v)

    #kernel values
    kerneldata = kernelvals(igd.op,tgeo,bgeo)

    #values & grad/div/curl of local shape functions
    tval = igd.test_local_space(tgeo) 
    bval = igd.trial_local_space(bgeo)

    #jacobian
    j = jacobian(tgeo) * jacobian(bgeo)
    
    integrand(igd.op, kerneldata,tval,tgeo,bval,tgeo) * j
end
=#


function integrand(viop::VSIESingleLayer, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:3]
    fy = @SVector[bvals[i].value for i in 1:4]

    dgx = @SVector[tvals[i].divergence for i in 1:3]
    dfy = @SVector[bvals[i].divergence for i in 1:4]

    G = kerneldata.green

    Ty = kerneldata.tau

    α = viop.α
    β = viop.β

    @SMatrix[α * dot(gx[i],Ty*fy[j]) * G + β * (dgx[i] * Ty*dfy[j])*G for i in 1:3, j in 1:4]
end

function integrand(viop::VSIEBoundary, kerneldata, tvals, tgeo, bvals, bgeo)

    dgx = @SVector[tvals[i].divergence for i in 1:3]
    fy = @SVector[bvals[i].value for i in 1:1]

    G = kerneldata.green

    Ty = kerneldata.tau

    α = viop.α

    @SMatrix[α * Ty * dgx[i] * fy[j] * G for i in 1:3, j in 1:1]
end

#=
function integrand(viop::VSIESingleLayer2, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:3]
    fy = @SVector[bvals[i].value for i in 1:4]

    dgx = @SVector[tvals[i].divergence for i in 1:3]
    dfy = @SVector[bvals[i].divergence for i in 1:4]

    G = kerneldata.green

    Ty = kerneldata.tau

    α = viop.α
    β = viop.β

    @SMatrix[α * dot(gx[i],Ty*fy[j]) * G - β * (dgx[i] * Ty*dfy[j])*G for i in 1:3, j in 1:4]
end


function integrand(viop::VSIEBoundaryT, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:1]
    dfy = @SVector[bvals[i].divergence for i in 1:3]

    G = kerneldata.green

    Tx = kerneldata.tau

    α = viop.α

    @SMatrix[α * Tx * gx[i] * dfy[j] * G for i in 1:3, j in 1:1]
end
=#

function integrand(viop::VSIEDoubleLayer, kerneldata, tvals, tgeo, bvals, bgeo)
    gx = @SVector[tvals[i].value for i in 1:3]
    fy = @SVector[bvals[i].value for i in 1:4]

    gradG = kerneldata.gradgreen
    
    Ty = kerneldata.tau
    
    α = viop.α

    @SMatrix[α * dot(cross(Ty*fy[j],gx[i]),gradG) for i in 1:3, j in 1:4]
end

#=
function integrand(viop::VSIEDoubleLayerT, kerneldata, tvals, tgeo, bvals, bgeo)
    gx = @SVector[tvals[i].value for i in 1:4]
    fy = @SVector[bvals[i].value for i in 1:3]


    gradG = kerneldata.gradgreen
 
    Ty = kerneldata.tau
    
    α = viop.α

    @SMatrix[α * transpose(cross(fy[j],gx[i])) *gradG for i in 1:4, j in 1:3]
end
=#

defaultquadstrat(op::VSIEOperator, tfs, bfs) = SauterSchwab3DQStrat(3,3,3,3,3,3)


function quaddata(op::VSIEOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts, trial_charts, qs::SauterSchwab3DQStrat)

    #The combinations of rules (6,7) and (5,7 are) BAAAADDDD
    # they result in many near singularity evaluations with any
    # resemblence of accuracy going down the drain! Simply don't!
    # (same for (5,7) btw...).
    t_qp = quadpoints(test_local_space,  test_charts,  (qs.outer_rule))
    b_qp = quadpoints(trial_local_space, trial_charts, (qs.inner_rule))

   
    sing_qp = (SauterSchwab3D._legendre(qs.sauter_schwab_1D,0,1), 
               SauterSchwab3D._shunnham2D(qs.sauter_schwab_2D),
               SauterSchwab3D._shunnham3D(qs.sauter_schwab_3D),
               SauterSchwab3D._shunnham4D(qs.sauter_schwab_4D),)


    return (tpoints=t_qp, bpoints=b_qp, sing_qp=sing_qp)
end

quadrule(op::VolumeSurfaceOperator, g::RefSpace, f::RefSpace, i, τ, j, σ, qd, qs) = qr_volume(op, g, f, i, τ, j, σ, qd, qs)


function qr_volume(op::VolumeSurfaceOperator, g::RefSpace, f::RefSpace, i, τ, j, σ, qd,
    qs::SauterSchwab3DQStrat)
  
    @assert (length(τ.vertices)==3 && length(σ.vertices)==4)  "Expected simplex wrong"

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
            dmin2 = min(dmin2, d2)
            if d2 < dtol
                push!(idx_t,j)
                push!(idx_s,i)
                hits +=1
                break
            end
        end
    end

    #singData = SauterSchwab3D.Singularity{D,hits}(idx_t, idx_s )
   
    hits == 3 && return SauterSchwab3D.CommonFace5D_S(SauterSchwab3D.Singularity5DFace(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
    hits == 2 && return SauterSchwab3D.CommonEdge5D_S(SauterSchwab3D.Singularity5DEdge(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
    hits == 1 && return SauterSchwab3D.CommonVertex5D_S(SauterSchwab3D.Singularity5DPoint(idx_t,idx_s),(qd.sing_qp[3],qd.sing_qp[2]))
    #hits == 0 && return SauterSchwab3D.PositiveDistance5D_S(SauterSchwab3D.Singularity5DPositiveDistance(),(qd.sing_qp[3],qd.sing_qp[2]))
    
    return DoubleQuadRule(
        qd[1][1,i],
        qd[2][1,j])

end

quadrule(op::BoundarySurfaceOperator, g::RefSpace, f::RefSpace, i, τ, j, σ, qd, qs) = qr_boundary(op, g, f, i, τ, j, σ, qd, qs)

function qr_boundary(op::BoundarySurfaceOperator, g::RefSpace, f::RefSpace, i, τ, j,  σ, qd,
    qs::SauterSchwab3DQStrat)
    
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
            dmin2 = min(dmin2, d2)
            if d2 < dtol
                push!(idx_t,i)
                push!(idx_s,j)
                hits +=1
                break
            end
        end
    end

    #singData = SauterSchwab3D.Singularity{D,hits}(idx_t, idx_s )
   

    hits == 3 && return SauterSchwab3D.CommonFace4D_S(SauterSchwab3D.Singularity4DFace(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[3]))
    hits == 2 && return SauterSchwab3D.CommonEdge4D_S(SauterSchwab3D.Singularity4DEdge(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2]))
    hits == 1 && return SauterSchwab3D.CommonVertex4D_S(SauterSchwab3D.Singularity4DPoint(idx_t,idx_s),(qd.sing_qp[2]))


    return DoubleQuadRule(
        qd[1][1,i],
        qd[2][1,j])

end

