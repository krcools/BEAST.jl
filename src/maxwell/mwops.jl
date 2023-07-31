abstract type MaxwellOperator3D <: IntegralOperator end
abstract type MaxwellOperator3DReg <: MaxwellOperator3D end

struct KernelValsMaxwell3D{T,U,P,Q}
    "gamma = im * wavenumber"
    gamma::U
    vect::P
    dist::T
    green::U
    gradgreen::Q
end

const inv_4pi = 1/(4pi)
function kernelvals(biop::MaxwellOperator3D, p, q)

    γ = biop.gamma
    r = cartesian(p) - cartesian(q)
    T = eltype(r)
    R = norm(r)
    γR = γ*R

    inv_R = 1/R

    expn = exp(-γR)
    green = expn * inv_R * T(inv_4pi)
    gradgreen = -(γ + inv_R) * green * inv_R * r

    KernelValsMaxwell3D(γ, r, R, green, gradgreen)
end


struct MWSingleLayer3D{T,U} <: MaxwellOperator3D
  gamma::T
  α::U
  β::U
end

MWSingleLayer3D(gamma)  = MWSingleLayer3D(gamma, -gamma, -1/(gamma))
MWWeaklySingular(gamma) = MWSingleLayer3D(gamma, 1, 0)
MWHyperSingular(gamma)  = MWSingleLayer3D(gamma, 0, 1)



export Maxwell3D

struct MWSingleLayer3DReg{T,U} <: MaxwellOperator3DReg
    gamma::T
    α::U
    β::U
end

struct MWSingleLayer3DSng{T,U} <: MaxwellOperator3D
    gamma::T
    α::U
    β::U
end

scalartype(op::MaxwellOperator3D) = typeof(op.gamma)

regularpart(op::MWSingleLayer3D) = MWSingleLayer3DReg(op.gamma, op.α, op.β)
singularpart(op::MWSingleLayer3D) = MWSingleLayer3DSng(op.gamma, op.α, op.β)


function _legendre(n,a,b)
    x, w = FastGaussQuadrature.gausslegendre(n)
    w .*= (b-a)/2
    x .= (x.+1)/2*(b-a).+a
    collect(zip(x,w))
end

defaultquadstrat(op::MaxwellOperator3D, tfs::Space, bfs::Space) = DoubleNumWiltonSauterQStrat(2,3,6,7,5,5,4,3)
defaultquadstrat(op::MaxwellOperator3D, tfs::RefSpace, bfs::RefSpace) = DoubleNumWiltonSauterQStrat(2,3,6,7,5,5,4,3)


function quaddata(op::MaxwellOperator3D,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts, trial_charts, qs::DoubleNumWiltonSauterQStrat)

    T = coordtype(test_charts[1])

    tqd = quadpoints(test_local_space,  test_charts,  (qs.outer_rule_far,qs.outer_rule_near))
    bqd = quadpoints(trial_local_space, trial_charts, (qs.inner_rule_far,qs.inner_rule_near))
     
    leg = (
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_vert,0,1)),
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_edge,0,1)),
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_face,0,1)),)

    return (tpoints=tqd, bpoints=bqd, gausslegendre=leg)
end

struct MWDoubleLayer3D{T} <: MaxwellOperator3D
  gamma::T
end

struct MWDoubleLayer3DSng{T} <: MaxwellOperator3D
    gamma::T
end

struct MWDoubleLayer3DReg{T} <: MaxwellOperator3DReg
    gamma::T
end

regularpart(op::MWDoubleLayer3D) = MWDoubleLayer3DReg(op.gamma)
singularpart(op::MWDoubleLayer3D) = MWDoubleLayer3DSng(op.gamma)

const LinearRefSpaceTriangle = Union{RTRefSpace, NDRefSpace, BDMRefSpace, NCrossBDMRefSpace}

function quadrule(op::MaxwellOperator3D, g::LinearRefSpaceTriangle, f::LinearRefSpaceTriangle,  i, τ, j, σ, qd,
      qs::DoubleNumWiltonSauterQStrat)

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

    h2 = volume(σ)
    xtol2 = 0.2 * 0.2
    k2 = abs2(op.gamma)
    max(dmin2*k2, dmin2/16h2) < xtol2 && return WiltonSERule(
        qd.tpoints[2,i],
        DoubleQuadRule(
            qd.tpoints[2,i],
            qd.bpoints[2,j],),)
    return DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j],)
end

function quadrule(op::MaxwellOperator3D, g::NDRefSpace, f::RTRefSpace,  i, τ, j, σ, qd,
  qs::DoubleNumWiltonSauterQStrat)

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

h2 = volume(σ)
xtol2 = 0.2 * 0.2
k2 = abs2(op.gamma)
max(dmin2*k2, dmin2/16h2) < xtol2 && return WiltonSERule(
    qd.tpoints[2,i],
    DoubleQuadRule(
        qd.tpoints[2,i],
        qd.bpoints[2,j],),)
return DoubleQuadRule(
    qd.tpoints[1,i],
    qd.bpoints[1,j],)
end

function quadrule(op::MaxwellOperator3D, g::BDMRefSpace, f::BDMRefSpace,  i, τ, j, σ, qd,
  qs::DoubleNumWiltonSauterQStrat)

  hits = 0
  dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))
  dmin2 = floatmax(eltype(eltype(τ.vertices)))
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

  h2 = volume(σ)
  xtol2 = 0.2 * 0.2
  k2 = abs2(op.gamma)
  return DoubleQuadRule(
      qd.tpoints[1,i],
      qd.bpoints[1,j],)
end



function qrib(op::MaxwellOperator3D, g::RTRefSpace, f::RTRefSpace, i, τ, j, σ, qd)

  dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))
  xtol = 0.2

  k = norm(op.gamma)

  hits = 0
  xmin = xtol
  for t in τ.vertices
    for s in σ.vertices
      d = norm(t-s)
      xmin = min(xmin, k*d)
      if d < dtol
        hits +=1
        break
      end
    end
  end

  hits == 3   && return BogaertSelfPatchStrategy(5)
  hits == 2   && return BogaertEdgePatchStrategy(8, 4)
  hits == 1   && return BogaertPointPatchStrategy(2, 3)
  rmin = xmin/k
  xmin < xtol && return WiltonSERule(
    qd.tpoints[1,i],
    DoubleQuadRule(
      qd.tpoints[2,i],
      qd.bpoints[2,j],
    ),
  )
  return DoubleQuadRule(
    qd.tpoints[1,i],
    qd.bpoints[1,j],
  )

end


function qrdf(op::MaxwellOperator3D, g::RTRefSpace, f::RTRefSpace, i, τ, j, σ, qd)
  # defines coincidence of points
  dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))

  # decides on whether to use singularity extraction
  xtol = 0.2

  k = norm(op.gamma)

  hits = 0
  xmin = xtol
  for t in τ.vertices
    for s in σ.vertices
      d = norm(t-s)
      xmin = min(xmin, k*d)
      if d < dtol
        hits +=1
        break
      end
    end
  end

  xmin < xtol && return WiltonSERule(
    qd.tpoints[1,i],
    DoubleQuadRule(
      qd.tpoints[2,i],
      qd.bpoints[2,j],
    ),
  )
  return DoubleQuadRule(
    qd.tpoints[1,i],
    qd.bpoints[1,j],
  )

end

################################################################################
#
#  Kernel definitions
#
################################################################################

const i4pi = 1 / (4pi)
function (igd::Integrand{<:MWSingleLayer3D})(x,y,f,g)
    α = igd.operator.α
    β = igd.operator.β
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1 / R
    green = exp(-γ*R)*(i4pi*iR)

    αG = α * green
    βG = β * green

    _integrands(f,g) do fi,gj
        αG * dot(fi.value, gj.value) + βG * dot(fi.divergence, gj.divergence)
    end
end

function (igd::Integrand{<:MWSingleLayer3DReg})(x,y,f,g)
    α = igd.operator.α
    β = igd.operator.β
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    γR = γ*R
    # iR = 1 / R
    green = (expm1(-γR) + γR - 0.5*γR^2) / (4pi*R)

    αG = α * green
    βG = β * green

    _integrands(f,g) do fi,gj
        αG * dot(fi.value, gj.value) + βG * dot(fi.divergence, gj.divergence)
    end
end


function (igd::Integrand{<:MWDoubleLayer3D})(x,y,f,g)
    
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)
    gradgreen = -(γ + iR) * green * (iR * r)

    fvalue = getvalue(f)
    gvalue = getvalue(g)
    G = cross.(Ref(gradgreen), gvalue)
    return _krondot(fvalue, G)
end


function (igd::Integrand{<:MWDoubleLayer3DReg})(x,y,f,g)
    
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    γR = γ*R
    iR = 1/R
    expo = exp(-γR)
    green = (expo - 1 + γR - 0.5*γR^2) * (i4pi*iR)
    gradgreen = ( -(γR + 1)*expo + (1 - 0.5*γR^2) ) * (i4pi*iR^3) * r

    fvalue = getvalue(f)
    gvalue = getvalue(g)
    G = cross.(Ref(gradgreen), gvalue)
    return _krondot(fvalue, G)
end


