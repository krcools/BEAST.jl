abstract type MaxwellOperator3D{T,K} <: IntegralOperator end
abstract type MaxwellOperator3DReg{T,K} <: MaxwellOperator3D{T,K} end

import base: ×, ⋅

scalartype(op::MaxwellOperator3D{T,K}) where {T, K <: Nothing} = T
scalartype(op::MaxwellOperator3D{T,K}) where {T, K} = promote_type(T, K)

gamma(op::MaxwellOperator3D{T,K}) where {T, K <: Nothing} = T(0)
gamma(op::MaxwellOperator3D{T,K}) where {T, K} = op.gamma

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

    γ = gamma(biop)
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
struct MWDoubleLayer3D{T,K} <: MaxwellOperator3D{T,K}
  alpha::T
  gamma::K
end

struct MWSingleLayer3D{T,U} <: MaxwellOperator3D{T,U}
  gamma::T
  α::U
  β::U
end

# struct TdGdn_ydj{T,U} <: MaxwellOperator3D{T,U}
#   gamma::T
#   α::U
# end

struct MWgreenint{T,U} <: MaxwellOperator3D{T,U}#check
  gamma::T
  α::U
end
struct MWgradgreendot{T,U} <: MaxwellOperator3D{T,U}#check
  gamma::T
  α::U
end

struct nXMWgreenint{T,U} <: MaxwellOperator3D{T,U}#check
  gamma::T
  α::U
end

struct nXdoublelayer{T,U} <: MaxwellOperator3D{T,U}#check
  gamma::T
  α::U
end

struct MWgreenintdotn{T,U}<: MaxwellOperator3D{T,U}#check
  gamma::T
  α::U
end
struct nXMWgreenintdotn{T,U}<: MaxwellOperator3D{T,U}#check
  gamma::T
  α::U
end
struct nXMWgradgreendot{T,U}<: MaxwellOperator3D{T,U}#check
  gamma::T
  α::U
end
struct MWgradgreendotn{T,U}<: MaxwellOperator3D{T,U}#check
  gamma::T
  α::U
end
struct MWdoublelayerXn{T,U}<: MaxwellOperator3D{T,U}#check
  gamma::T
  α::U
end
struct nXMWdoublelayerXn{T,U}<: MaxwellOperator3D{T,U}#check haakjes rechts
  gamma::T
  α::U
end
struct ndotMWdoublelayer{T,U}<: MaxwellOperator3D{T,U}#check
  gamma::T
  α::U
end
struct ndotMWgreenint{T,U}<: MaxwellOperator3D{T,U}#cehck
  gamma::T
  α::U
end
struct ndotMWgreenintdotn{T,U}<: MaxwellOperator3D{T,U}#check
  gamma::T
  α::U
end
struct ndotMWgradgreen{T,U}<: MaxwellOperator3D{T,U}#check
  gamma::T
  α::U
end

×(n::NormalVector,op::MWgreenint) = nXMWgreenint(op.gamma,op.α)
×(n::NormalVector,op::MWDoubleLayer3D) = nXdoublelayer(op.gamma,op.α)
⋅(op::MWgreenint,n::NormalVector) = MWgreenintdotn(op.gamma,op.α)
⋅(op::nXMWgreenint,n::NormalVector) = nXMWgreenintdotn(op.gamma,op.α)
×(n::NormalVector,op::MWgreenintdotn) = nXMWgreenintdotn(op.gamma,op.α)
×(n::NormalVector,op::MWgradgreendot) = nXMWgradgreendot(op.gamma,op.α)
⋅(op::MWgradgreendot,n::NormalVector) = MWgradgreendotn(op.gamma,op.α)
×(op::MWDoubleLayer3D,n::NormalVector) = MWdoublelayerXn(op.gamma,op.α)
×(n::NormalVector,op::MWdoublelayerXn) = nXMWdoublelayerXn(op.gamma,op.α)
⋅(n::NormalVector,op::MWDoubleLayer3D) = ndotMWdoublelayer(op.gamma,op.α)
⋅(n::NormalVector,op::MWgreenint) = ndotMWgreenint(op.gamma,op.α)
⋅(n::NormalVector,op::MWgreenintdotn) = ndotMWgreenintdotn(op.gamma,op.α)
⋅(op::ndotMWgreenint,n::NormalVector)= ndotMWgreenintdotn(op.gamma,op.α)
⋅(n::NormalVector,op::MWgradgreen) = ndotMWgradgreen(op.gamma,op.α)

# struct MWgradgreenint{T,U} <: MaxwellOperator3D{T,U}
#   gamma::T
#   α::U
# end
# struct MWngradgreenint{T,U} <: MaxwellOperator3D{T,U}
#   gamma::T
#   α::U
# end
scalartype(op::MWSingleLayer3D{T,U}) where {T,U} = promote_type(T,U)
sign_upon_permutation(op::MWSingleLayer3D, I, J) = 1

MWSingleLayer3D(gamma)  = MWSingleLayer3D(gamma, -gamma, -1/(gamma))
MWWeaklySingular(gamma) = MWSingleLayer3D(gamma, 1, 0)
MWHyperSingular(gamma)  = MWSingleLayer3D(gamma, 0, 1)



export Maxwell3D

struct MWSingleLayer3DReg{T,U} <: MaxwellOperator3DReg{T,U}
    gamma::T
    α::U
    β::U
end

struct MWSingleLayer3DSng{T,U} <: MaxwellOperator3D{T,U}
    gamma::T
    α::U
    β::U
end

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



# struct MWnDoubleLayer3D{T,K} <: MaxwellOperator3D{T,K}
#   alpha::T
#   gamma::K
# end

sign_upon_permutation(op::MWDoubleLayer3D, I, J) = 1

struct MWDoubleLayer3DSng{T,K} <: MaxwellOperator3D{T,K}
    alpha::T
    gamma::K
end

struct MWDoubleLayer3DReg{T,K} <: MaxwellOperator3DReg{T,K}
    alpha::T
    gamma::K
end

MWDoubleLayer3D(gamma) = MWDoubleLayer3D(1.0, gamma) # For legacy purposes
MWnDoubleLayer3D(gamma) = MWnDoubleLayer3D(1.0, gamma) # For legacy purposes

regularpart(op::MWDoubleLayer3D) = MWDoubleLayer3DReg(op.alpha, op.gamma)
singularpart(op::MWDoubleLayer3D) = MWDoubleLayer3DSng(op.alpha, op.gamma)

function quadrule(op::MaxwellOperator3D, g::RTRefSpace, f::RTRefSpace,  i, τ, j, σ, qd,
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
    k2 = abs2(gamma(op))
    println(max(dmin2*k2, dmin2/16h2) < xtol2)
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
  k2 = abs2(gamma(op))
  return DoubleQuadRule(
      qd.tpoints[1,i],
      qd.bpoints[1,j],)
end


function qrib(op::MaxwellOperator3D, g::RTRefSpace, f::RTRefSpace, i, τ, j, σ, qd)

  dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))
  xtol = 0.2

  k = norm(gamma(op))

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

  k = norm(gamma(op))

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
ttrace!(op::MWSingleLayer3D, mesh, or) = ZeroOperator(), SameBase()
strace!(op::MWSingleLayer3D, mesh, or) = ZeroOperator(), mesh

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
function (igd::Integrand{<:ndotMWdoublelayer})(x,y,f,g)
    
  γ = igd.operator.gamma

  r = cartesian(x) - cartesian(y)
  R = norm(r)
  iR = 1/R
  green = exp(-γ*R)*(iR*i4pi)
  gradgreen = -(γ + iR) * green * (iR * r)
  nx = normal(x)
  fvalue = getvalue(f)
  gvalue = getvalue(g)
  G = dot.(Ref(nx),cross.(Ref(gradgreen), gvalue))
  return _krondot(fvalue, G)
end
function (igd::Integrand{<:MWdoublelayerXn})(x,y,f,g)
    
  γ = igd.operator.gamma

  r = cartesian(x) - cartesian(y)
  R = norm(r)
  iR = 1/R
  green = exp(-γ*R)*(iR*i4pi)
  gradgreen = -(γ + iR) * green * (iR * r)
  ny = normal(y)
  fvalue = getvalue(f)
  gvalue = getvalue(g)
  G = cross.(Ref(cross(gradgreen,n)), gvalue)
  return _krondot(fvalue, G)
end

function (igd::Integrand{<:nXdoublelayer})(x,y,f,g)
    
  γ = igd.operator.gamma

  r = cartesian(x) - cartesian(y)
  R = norm(r)
  iR = 1/R
  green = exp(-γ*R)*(iR*i4pi)
  gradgreen = -(γ + iR) * green * (iR * r)
  nx = normal(x)
  fvalue = getvalue(f)
  gvalue = getvalue(g)
  G = cross.(Ref(cross(nx,gradgreen)), gvalue)
  return _krondot(fvalue, G)
end
function (igd::Integrand{<:nXMWdoublelayerXn})(x,y,f,g)
    
  γ = igd.operator.gamma

  r = cartesian(x) - cartesian(y)
  R = norm(r)
  iR = 1/R
  green = exp(-γ*R)*(iR*i4pi)
  gradgreen = -(γ + iR) * green * (iR * r)
  nx = normal(x)
  ny = normal(y)
  fvalue = getvalue(f)
  gvalue = getvalue(g)
  G = Ref(cross(nx,cross(gradgreen,ny))).* gvalue
  return _krondot(fvalue, G)
end
# ttrace!(op::MWDoubleLayer3D, mesh, orientation::Inside) = 1/2*NCross(), SameBase()
# strace!(op::MWDoubleLayer3D, mesh, orientation::Inside) = 1/2*NCross(), mesh
# ttrace!(op::MWDoubleLayer3D, mesh, orientation::Outside) = -1/2*NCross(), SameBase()
# strace!(op::MWDoubleLayer3D, mesh, orientation::Outside) = -1/2*NCross(), mesh

# function (igd::Integrand{<:MWnDoubleLayer3D})(x,y,f,g)
    
#   γ = igd.operator.gamma

#   r = cartesian(x) - cartesian(y)
#   R = norm(r)
#   iR = 1/R
#   green = exp(-γ*R)*(iR*i4pi)
#   gradgreen = -(γ + iR) * green * (iR * r)
#   normaly = normal(y)
#   fvalue = getvalue(f)
#   gvalue = getvalue(g)
#   t1 = Ref(normaly).*gvalue
#   G = cross.(Ref(gradgreen), t1)
#   return _krondot(fvalue, G)
# end

# ttrace!(op::MWnDoubleLayer3D, mesh, or) = ZeroOperator(), SameBase()
# strace!(op::MWnDoubleLayer3D, mesh, or) = ZeroOperator(), mesh

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

# function (igd::Integrand{<:MWngreenint})(x,y,f,g)
#     γ = igd.operator.gamma

#     r = cartesian(x) - cartesian(y)
#     R = norm(r)
#     γR = γ*R
#     iR = 1/R
#     expo = exp(-γR)
#     green = (expo - 1 + γR - 0.5*γR^2) * (i4pi*iR)
    
#     normaly = normal(y)
#     fvalue = getvalue(f)
#     gvalue = getvalue(g)
#     t = normaly*green
#     t2 = Ref(t).*gvalue
#     return _krondot(fvalue,t2)
# end

# ttrace!(op::MWngreenint, mesh, or) = ZeroOperator(), SameBase()
# strace!(op::MWngreenint, mesh, or) = ZeroOperator(), mesh



function (igd::Integrand{<:MWgreenint})(x,y,f,g)
  γ = igd.operator.gamma

  r = cartesian(x) - cartesian(y)
  R = norm(r)
  iR = 1/R
  green = exp(-γ*R)*(iR*i4pi)
  gradgreen = -(γ + iR) * green * (iR * r)

  fvalue = getvalue(f)
  gvalue = getvalue(g)

  t2 = green*gvalue
  return _krondot(fvalue,t2)
end
function (igd::Integrand{<:ndotMWgreenint})(x,y,f,g)
  γ = igd.operator.gamma

  r = cartesian(x) - cartesian(y)
  R = norm(r)
  iR = 1/R
  green = exp(-γ*R)*(iR*i4pi)
  gradgreen = -(γ + iR) * green * (iR * r)
  nx = normal(x)
  fvalue = getvalue(f)
  gvalue = getvalue(g)

  t2 = dot.(Ref(nx),green*gvalue)
  return _krondot(fvalue,t2)
end
function (igd::Integrand{<:ndotMWgreenintdotn})(x,y,f,g)
  γ = igd.operator.gamma

  r = cartesian(x) - cartesian(y)
  R = norm(r)
  iR = 1/R
  green = exp(-γ*R)*(iR*i4pi)
  gradgreen = -(γ + iR) * green * (iR * r)
  nx = normal(x)
  ny = normal(y)
  fvalue = getvalue(f)
  gvalue = getvalue(g)

  t2 = dot.(Ref(nx),Ref(ny))*green*gvalue
  return _krondot(fvalue,t2)
end
function (igd::Integrand{<:MWgreenintdotn})(x,y,f,g)
  γ = igd.operator.gamma

  r = cartesian(x) - cartesian(y)
  R = norm(r)
  iR = 1/R
  green = exp(-γ*R)*(iR*i4pi)
  gradgreen = -(γ + iR) * green * (iR * r)
  ny = normal(y)
  fvalue = getvalue(f)
  gvalue = getvalue(g)

  t2 = green*gvalue.*Ref(ny)
  return _krondot(fvalue,t2)
end
function (igd::Integrand{<:nXMWgreenintdotn})(x,y,f,g)
  γ = igd.operator.gamma

  r = cartesian(x) - cartesian(y)
  R = norm(r)
  iR = 1/R
  green = exp(-γ*R)*(iR*i4pi)
  gradgreen = -(γ + iR) * green * (iR * r)
  ny = normal(y)
  nx = normal(x)
  fvalue = getvalue(f)
  gvalue = getvalue(g)
  t1 = green*gvalue.*Ref(ny)
  t2 = cross.(Ref(nx),t1)
  return _krondot(fvalue,t2)
end

function (igd::Integrand{<:nXMWgreenint})(x,y,f,g)
  γ = igd.operator.gamma

  r = cartesian(x) - cartesian(y)
  R = norm(r)
  iR = 1/R
  green = exp(-γ*R)*(iR*i4pi)
  gradgreen = -(γ + iR) * green * (iR * r)
  nx = normal(x)
  fvalue = getvalue(f)
  gvalue = getvalue(g)

  t2 = cross.(Ref(nx),green*gvalue)
  return _krondot(fvalue,t2)
end

function (igd::Integrand{<:MWgradgreendot})(x,y,f,g)
  γ = igd.operator.gamma

  r = cartesian(x) - cartesian(y)
  R = norm(r)
  iR = 1/R
  green = exp(-γ*R)*(iR*i4pi)
  gradgreen = -(γ + iR) * green * (iR * r)

  fvalue = getvalue(f)
  gvalue = getvalue(g)

  t2 = dot.(gvalue,Ref(gradgreen))
  return _krondot(fvalue,t2)
end
function (igd::Integrand{<:ndotMWgradgreen})(x,y,f,g)
  γ = igd.operator.gamma

  r = cartesian(x) - cartesian(y)
  R = norm(r)
  iR = 1/R
  green = exp(-γ*R)*(iR*i4pi)
  gradgreen = -(γ + iR) * green * (iR * r)
  nx = normal(x)
  fvalue = getvalue(f)
  gvalue = getvalue(g)

  t2 = gvalue*dot(nx,gradgreen)
  return _krondot(fvalue,t2)
end

function (igd::Integrand{<:MWgradgreendotn})(x,y,f,g)
  γ = igd.operator.gamma

  r = cartesian(x) - cartesian(y)
  R = norm(r)
  iR = 1/R
  green = exp(-γ*R)*(iR*i4pi)
  gradgreen = -(γ + iR) * green * (iR * r)

  fvalue = getvalue(f)
  gvalue = getvalue(g)
  ny = normal(y)
  t2 = gvalue * dot(ny,Ref(gradgreen))
  return _krondot(fvalue,t2)
end
function (igd::Integrand{<:nXMWgradgreendot})(x,y,f,g)
  γ = igd.operator.gamma

  r = cartesian(x) - cartesian(y)
  R = norm(r)
  iR = 1/R
  green = exp(-γ*R)*(iR*i4pi)
  gradgreen = -(γ + iR) * green * (iR * r)
  nx = normal(x)
  fvalue = getvalue(f)
  gvalue = getvalue(g)

  t2 = gvalue.*Ref(cross(nx,gradgreen))
  return _krondot(fvalue,t2)
end

# ttrace!(op::MWgreenint, mesh, or) = ZeroOperator(), SameBase()
# strace!(op::MWgreenint, mesh, or) = ZeroOperator(), mesh
# ntrace!(op::MWgreenint, mesh, or) = ZeroOperator(), mesh

# function (igd::Integrand{<:MWgradgreenint})(x,y,f,g)
#   γ = igd.operator.gamma

#   r = cartesian(x) - cartesian(y)
#   R = norm(r)
#   γR = γ*R
#   iR = 1/R
#   expo = exp(-γR)
#   green = (expo - 1 + γR - 0.5*γR^2) * (i4pi*iR)
#   gradgreen = ( -(γR + 1)*expo + (1 - 0.5*γR^2) ) * (i4pi*iR^3) * r
  
#   fvalue = getvalue(f)
#   gvalue = getvalue(g)

#   t2 = dot.(gvalue,Ref(gradgreen))

#   return _krondot(fvalue,t2)
# end
# ttrace!(op::MWgradgreenint, mesh, or) = ZeroOperator(), SameBase()
# strace!(op::MWgradgreenint, mesh, or) = ZeroOperator(), mesh
# ntrace!(op::MWgradgreenint, mesh, or) = 1/2*Identity(), mesh

# function (igd::Integrand{<:MWngradgreenint})(x,y,f,g)
#   γ = igd.operator.gamma

#   r = cartesian(x) - cartesian(y)
#   R = norm(r)
#   γR = γ*R
#   iR = 1/R
#   expo = exp(-γR)
#   green = (expo - 1 + γR - 0.5*γR^2) * (i4pi*iR)
#   gradgreen = ( -(γR + 1)*expo + (1 - 0.5*γR^2) ) * (i4pi*iR^3) * r
  
#   fvalue = getvalue(f)
#   gvalue = getvalue(g)
#   normaly = normal(y)
#   t1 = Ref(normaly).*gvalue
#   t2 = dot.(t1,Ref(gradgreen))

#   return _krondot(fvalue,t2)
# end

# scallartrace!(op::MWngradgreenint, mesh, or) = 1/2*Identity(), SameBase()
################################################################################
#
#  Handling of operator parameters (Helmholtz and Maxwell)
#
################################################################################

function gamma_wavenumber_handler(gamma, wavenumber)
  if (gamma !== nothing) && (wavenumber !== nothing)
      error("Supply one of (not both) gamma or wavenumber")
  end

  if gamma === nothing && (wavenumber !== nothing)
      if iszero(real(wavenumber))
          gamma = -imag(wavenumber)
      else
          gamma = im*wavenumber
      end
  end

  return gamma, wavenumber
end

function operator_parameter_handler(alpha, gamma, wavenumber)

  gamma, wavenumber = gamma_wavenumber_handler(gamma, wavenumber)

  if alpha === nothing
      if gamma !== nothing
          alpha = one(real(typeof(gamma)))
      else
          # We are dealing with a static problem. Default to double precision.
          alpha = 1.0
      end
  end

  return alpha, gamma
end