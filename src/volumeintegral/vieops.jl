abstract type VIEOperator <: IntegralOperator end
abstract type VolumeOperator <: VIEOperator end
abstract type BoundaryOperator <: VIEOperator end

struct KernelValsVIE{T,U,P,Q,K}
    gamma::U
    vect::P
    dist::T
    green::U
    gradgreen::Q
    tau::K
end

function kernelvals(viop::VIEOperator, p ,q)
    Y = viop.gamma;
    r = cartesian(p)-cartesian(q)
    R = norm(r)
    yR = Y*R

    expn = exp(-yR)
    green = expn / (4*pi*R)
    gradgreen = - (Y +1/R)*green/R*r  # Derivation after p (test variable)
    

    tau = viop.tau(cartesian(q))

    KernelValsVIE(Y,r,R, green, gradgreen,tau)
end

struct VIESingleLayer{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    β::U
    tau::P
end

struct VIEBoundary{T,U,P} <: BoundaryOperator
    gamma::T
    α::U
    tau::P
end

struct VIESingleLayer2{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    β::U
    tau::P
end

struct VIEBoundary2{T,U,P} <: BoundaryOperator
    gamma::T
    α::U
    tau::P
end


struct VIEDoubleLayer{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    tau::P
end


struct VIEhhVolume{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    tau::P
end

struct VIEhhBoundary{T,U,P} <: BoundaryOperator
    gamma::T
    α::U
    tau::P
end

struct VIEhhVolumegradG{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    tau::P
end

struct VIEhhVolumek0{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    tau::P
end


scalartype(op::VIEOperator) = typeof(op.gamma)

export VIE

struct VIEIntegrand{S,T,O,K,L}
    test_tetrahedron_element::S
    trial_tetrahedron_element::T
    op::O
    test_local_space::K
    trial_local_space::L
end


function (igd::VIEIntegrand)(u,v)

    #mesh points
    tgeo = neighborhood(igd.test_tetrahedron_element,v)
    bgeo = neighborhood(igd.trial_tetrahedron_element,u)

    #kernel values
    kerneldata = kernelvals(igd.op,tgeo,bgeo)

    #values & grad/div/curl of local shape functions
    tval = igd.test_local_space(tgeo) 
    bval = igd.trial_local_space(bgeo)

    #jacobian
    j = jacobian(tgeo) * jacobian(bgeo)
    
    integrand(igd.op, kerneldata,tval,tgeo,bval,bgeo) * j
end


function integrand(viop::VIESingleLayer, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:4]
    fy = @SVector[bvals[i].value for i in 1:4]

    dgx = @SVector[tvals[i].divergence for i in 1:4]
    dfy = @SVector[bvals[i].divergence for i in 1:4]

    G = kerneldata.green
    gradG = kerneldata.gradgreen

    Ty = kerneldata.tau

    α = viop.α
    β = viop.β

    @SMatrix[α * dot(gx[i],Ty*fy[j]) * G - β * dot(dgx[i] * (Ty * fy[j]), gradG) for i in 1:4, j in 1:4]
end


function integrand(viop::VIESingleLayer2, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:6]
    fy = @SVector[bvals[i].value for i in 1:6]

    dgx = @SVector[tvals[i].curl for i in 1:6]
    dfy = @SVector[bvals[i].curl for i in 1:6]

    G = kerneldata.green
    gradG = kerneldata.gradgreen

    Ty = kerneldata.tau

    α = viop.α
    β = viop.β

 
    @SMatrix[dot(dgx[i] , cross(gradG, Ty*fy[j])) for i in 1:6, j in 1:6]
end

function integrand(viop::VIEBoundary, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:1]
    fy = @SVector[bvals[i].value for i in 1:4]

    gradG = kerneldata.gradgreen

    Ty = kerneldata.tau

    α = viop.α

    @SMatrix[α * gx[i] * dot(Ty * fy[j], gradG) for i in 1:1, j in 1:4]
end

function integrand(viop::VIEBoundary2, kerneldata, tvals, tgeo, bvals, bgeo)


    gx = @SVector[tvals[i].value for i in 1:3]
    fy = @SVector[bvals[i].value for i in 1:6]


    gradG = kerneldata.gradgreen

    Ty = kerneldata.tau

    α = viop.α

    @SMatrix[α * dot(gx[i] , cross(gradG, Ty*fy[j])) for i in 1:3, j in 1:6]
end

function integrand(viop::VIEDoubleLayer, kerneldata, tvals, tgeo, bvals, bgeo)
    gx = tvals[1]
    fy = bvals[1]

    gradG = kerneldata.gradgreen
    
    Ty = kerneldata.tau
    
    α = viop.α

    t = α * dot(gx, cross(gradG, Ty*fy))
end


# Integrands for the operators of the Lippmann Schwinger Volume Integral Equation:

function integrand(viop::VIEhhVolume, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:4]
    fy = @SVector[bvals[i].value for i in 1:4]

    dgx = @SVector[tvals[i].gradient for i in 1:4]
    dfy = @SVector[bvals[i].gradient for i in 1:4]

    G = kerneldata.green
    gradG = kerneldata.gradgreen

    Ty = kerneldata.tau

    α = viop.α

    return @SMatrix[α * dot(dgx[i], G*Ty*dfy[j]) for i in 1:4, j in 1:4]
end


function integrand(viop::VIEhhBoundary, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:3]
    dfy = @SVector[bvals[i].gradient for i in 1:4]

    G = kerneldata.green
    gradG = kerneldata.gradgreen

    Ty = kerneldata.tau

    α = viop.α

    return @SMatrix[α * dot( tgeo.patch.normals[1]*gx[i],G*Ty*dfy[j]) for i in 1:3, j in 1:4]
end

function integrand(viop::VIEhhVolumek0, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:4]
    fy = @SVector[bvals[i].value for i in 1:4]

    dgx = @SVector[tvals[i].gradient for i in 1:4]
    dfy = @SVector[bvals[i].gradient for i in 1:4]

    G = kerneldata.green
    gradG = kerneldata.gradgreen

    Ty = kerneldata.tau

    α = viop.α

    return @SMatrix[α * gx[i]*G*Ty*fy[j] for i in 1:4, j in 1:4]
end

function integrand(viop::VIEhhVolumegradG, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:4]
    fy = @SVector[bvals[i].value for i in 1:4]

    dgx = @SVector[tvals[i].gradient for i in 1:4]
    dfy = @SVector[bvals[i].gradient for i in 1:4]

    G = kerneldata.green
    gradG = -kerneldata.gradgreen # "-" to get nabla'G(r,r')

    Ty = kerneldata.tau

    α = viop.α

    return @SMatrix[α * gx[i] * dot(gradG, Ty*dfy[j]) for i in 1:4, j in 1:4]
end




defaultquadstrat(op::VIEOperator, tfs, bfs) = SauterSchwab3DQStrat(3,3,3,3,3,3)


function quaddata(op::VIEOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts, trial_charts, qs::SauterSchwab3DQStrat)

    #The combinations of rules (6,7) and (5,7 are) BAAAADDDD
    # they result in many near singularity evaluations with any
    # resemblence of accuracy going down the drain! Simply don't!
    # (same for (5,7) btw...).
    t_qp = quadpoints(test_local_space,  test_charts,  (qs.outer_rule,))
    b_qp = quadpoints(trial_local_space, trial_charts, (qs.inner_rule,))

   
    sing_qp = (SauterSchwab3D._legendre(qs.sauter_schwab_1D,0,1), 
               SauterSchwab3D._shunnham2D(qs.sauter_schwab_2D),
               SauterSchwab3D._shunnham3D(qs.sauter_schwab_3D),
               SauterSchwab3D._shunnham4D(qs.sauter_schwab_4D),)


    return (tpoints=t_qp, bpoints=b_qp, sing_qp=sing_qp)
end

quadrule(op::VolumeOperator, g::RefSpace, f::RefSpace, i, τ, j, σ, qd, qs) = qr_volume(op, g, f, i, τ, j, σ, qd, qs)


function qr_volume(op::VolumeOperator, g::RefSpace, f::RefSpace, i, τ, j, σ, qd,
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

quadrule(op::BoundaryOperator, g::RefSpace, f::RefSpace, i, τ, j, σ, qd, qs) = qr_boundary(op, g, f, i, τ, j, σ, qd, qs)

function qr_boundary(op::BoundaryOperator, g::RefSpace, f::RefSpace, i, τ, j,  σ, qd,
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

