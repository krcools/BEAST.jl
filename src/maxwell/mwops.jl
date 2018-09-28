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

function kernelvals(biop::MaxwellOperator3D, p, q)

    γ = biop.gamma
    r = cartesian(p) - cartesian(q)
    R = norm(r)
    γR = γ*R

    expn = exp(-γR)
    green = expn / (4pi*R)
    gradgreen = -(γ + 1/R) * green / R * r

    KernelValsMaxwell3D(γ, r, R, green, gradgreen)
end

function kernelvals(kernel::MaxwellOperator3DReg, p, q)

    γ = kernel.gamma
    r = p.cart - q.cart
    R = norm(r)
    γR = γ*R

    Exp = exp(-γ*R)
    green = (Exp - 1 + γR - 0.5*γR^2) / (4pi*R)
    gradgreen = ( - (γR + 1)*Exp + (1 - 0.5*γR^2) ) * (r/R^3) / (4π)

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
    x = (x.+1)/2*(b-a).+a
    collect(zip(x,w))
end

function quaddata(op::MaxwellOperator3D,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts, trial_charts)

    a, b = 0.0, 1.0
    # CommonVertex, CommonEdge, CommonFace rules

    tqd = quadpoints(test_local_space, test_charts, (2,6))
    bqd = quadpoints(trial_local_space, trial_charts, (3,7))
    leg = (_legendre(3,a,b), _legendre(4,a,b), _legendre(5,a,b),)

    # High accuracy rules (use them e.g. in LF MFIE scenarios)
    # tqd = quadpoints(test_local_space, test_charts, (8,8))
    # bqd = quadpoints(trial_local_space, trial_charts, (8,9))
    # leg = (_legendre(8,a,b), _legendre(10,a,b), _legendre(5,a,b),)


    return (tpoints=tqd, bpoints=bqd, gausslegendre=leg)
end




# use Union type so this code can be shared between the operator
# and its regular part.
MWSL3DGen = Union{MWSingleLayer3D,MWSingleLayer3DReg}
function integrand(biop::MWSL3DGen, kerneldata, tvals, tgeo, bvals, bgeo)

  gx = tvals[1]
  fy = bvals[1]

  dgx = tvals[2]
  dfy = bvals[2]

  G = kerneldata.green
  γ = kerneldata.gamma

  α = biop.α
  β = biop.β

  t = (α * dot(gx, fy) + β * (dgx*dfy)) * G
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

const MWDL3DGen = Union{MWDoubleLayer3D,MWDoubleLayer3DReg}
function integrand(biop::MWDL3DGen, kerneldata, tvals, tgeo, bvals, bgeo)
    g = tvals[1]
    f = bvals[1]
    ∇G = kerneldata.gradgreen
    g ⋅ (∇G × f)
end


quadrule(op::MaxwellOperator3D, g::RTRefSpace, f::RTRefSpace, i, τ, j, σ, qd) = qrss(op, g, f, i, τ, j, σ, qd)

function qrss(op, g, f, i, τ, j, σ, qd)
    # defines coincidence of points
    dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))

    # decides on whether to use singularity extraction
    xtol = 0.2
    k = norm(op.gamma)

    hits = 0
    dmin = floatmax(eltype(eltype(τ.vertices)))
    for t in τ.vertices
        for s in σ.vertices
            d = norm(t-s)
            dmin = min(dmin, d)
            if d < dtol
                hits +=1
                break
            end
        end
    end

    hits == 3 && return SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
    hits == 2 && return SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
    hits == 1 && return SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])

    h = sqrt(volume(σ))
    max(dmin*k, dmin/4h) < xtol && return WiltonSEStrategy(
        qd.tpoints[2,i],
        DoubleQuadStrategy(
            qd.tpoints[2,i],
            qd.bpoints[2,j],),)
    return DoubleQuadStrategy(
        qd.tpoints[1,i],
        qd.bpoints[1,j],)
end


function qrib(op::MaxwellOperator3D, g::RTRefSpace, f::RTRefSpace, i, τ, j, σ, qd)
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

    hits == 3   && return BogaertSelfPatchStrategy(5)
    hits == 2   && return BogaertEdgePatchStrategy(8, 4)
    hits == 1   && return BogaertPointPatchStrategy(2, 3)
    rmin = xmin/k
    xmin < xtol && return WiltonSEStrategy(
      qd.tpoints[1,i],
      DoubleQuadStrategy(
        qd.tpoints[2,i],
        qd.bpoints[2,j],
      ),
    )
    return DoubleQuadStrategy(
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

  xmin < xtol && return WiltonSEStrategy(
    qd.tpoints[1,i],
    DoubleQuadStrategy(
      qd.tpoints[2,i],
      qd.bpoints[2,j],
    ),
  )
  return DoubleQuadStrategy(
    qd.tpoints[1,i],
    qd.bpoints[1,j],
  )

end
