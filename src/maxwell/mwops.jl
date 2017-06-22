import Base.cross
import WiltonInts84

export MWSingleLayer3D, MWHyperSingular, MWWeaklySingular
export MWDoubleLayer3D
export PlaneWaveMW
export TangTraceMW, CrossTraceMW

export curl

abstract type MaxwellOperator3D <: IntegralOperator end
abstract type MaxwellOperator3DReg <: MaxwellOperator3D end

immutable KernelValsMaxwell3D{T,U,P,Q}
    "gamma = im * wavenumber"
    gamma::U
    #vect::Pt{3,T}
    vect::P
    dist::T
    green::U
    #gradgreen::Pt{3,U}
    gradgreen::Q
end

function kernelvals(biop::MaxwellOperator3D, p, q)

    γ = biop.gamma
    r = p.cart - q.cart
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

immutable MWSingleLayer3D{T,U} <: MaxwellOperator3D
  #wavenumber::T
  gamma::T
  α::U
  β::U
end

MWSingleLayer3D(gamma)  = MWSingleLayer3D(gamma, -gamma, -1/(gamma))
MWWeaklySingular(gamma) = MWSingleLayer3D(gamma, 1, 0)
MWHyperSingular(gamma)  = MWSingleLayer3D(gamma, 0, 1)

immutable MWSingleLayer3DReg{T,U} <: MaxwellOperator3DReg
    gamma::T
    α::U
    β::U
end

immutable MWSingleLayer3DSng{T,U} <: MaxwellOperator3D
    gamma::T
    α::U
    β::U
end

scalartype(op::MaxwellOperator3D) = typeof(op.gamma)

regularpart(op::MWSingleLayer3D) = MWSingleLayer3DReg(op.gamma, op.α, op.β)
singularpart(op::MWSingleLayer3D) = MWSingleLayer3DSng(op.gamma, op.α, op.β)


function quaddata(op::MaxwellOperator3D, g::RefSpace, f::RefSpace, tels, bels)

    tqd = quadpoints(g, tels, (2,6))
    bqd = quadpoints(f, bels, (3,7))

    return QuadData(tqd, bqd)
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

  #t1 = α * G * dot(gx, fy)
  #t2 = β * G * dot(dgx, dfy)
  t = (α * dot(gx, fy) + β * (dgx*dfy)) * G
  #return t1 + t2
  return t
end

immutable MWDoubleLayer3D{T} <: MaxwellOperator3D
  gamma::T
end

immutable MWDoubleLayer3DSng{T} <: MaxwellOperator3D
    gamma::T
end

immutable MWDoubleLayer3DReg{T} <: MaxwellOperator3DReg
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

function innerintegrals!(op::MWSingleLayer3DSng, p, g, f, t, s, z, strat::WiltonSEStrategy, dx)

    γ = op.gamma
    x = cartesian(p)
    n = cross(s[1]-s[3],s[2]-s[3])
    n /= norm(n)
    ρ = x - ((x-s[1]) ⋅ n) * n

    scal, vec = WiltonInts84.wiltonints(s[1], s[2], s[3], x, Val{1})

    # \int \frac{1}{4 \pi R}
    ∫G = (scal[2] - γ*scal[3] + 0.5*γ^2*scal[4]) / (4π)
    # \int \frac{y}{4 \pi R}
    ∫Gy = SVector((
        (vec[2][1] + scal[2]*ρ[1] - γ*(vec[3][1]+scal[3]*ρ[1]) + 0.5*γ^2*(vec[4][1]+scal[4]*ρ[1]))/(4π),
        (vec[2][2] + scal[2]*ρ[2] - γ*(vec[3][2]+scal[3]*ρ[2]) + 0.5*γ^2*(vec[4][2]+scal[4]*ρ[2]))/(4π),
        (vec[2][3] + scal[2]*ρ[3] - γ*(vec[3][3]+scal[3]*ρ[3]) + 0.5*γ^2*(vec[4][3]+scal[4]*ρ[3]))/(4π),
    ))

    c₁ = op.α
    c₂ = op.β

    α = 1 / volume(t) / volume(s) / 4
    for i in 1 : numfunctions(g)
        a = t[i]
        g = x - a
        dg = 2

        for j in 1 : numfunctions(f)
            b = s[j]

            ∫Gf = SVector(∫Gy[1]-∫G*b[1], ∫Gy[2]-∫G*b[2], ∫Gy[3]-∫G*b[3])
            ∫Gdf = 2 * ∫G
            dg∫Gf = g[1]*∫Gf[1] + g[2]*∫Gf[2] + g[3]*∫Gf[3]
            z[i,j] += ( α*c₁*dg∫Gf + α*c₂*dg*∫Gdf ) * dx

        end # next j
    end #

end


function innerintegrals!(op::MWDoubleLayer3DSng, p, g, f, t, s, z, strat::WiltonSEStrategy, dx)

    γ = op.gamma
    x = cartesian(p)
    n = cross(s[1]-s[3],s[2]-s[3])
    n /= norm(n)
    ρ = x - ((x-s[1]) ⋅ n) * n

    #scal, vec, grad = wiltonints(s[1], s[2], s[3], x)
    scal, vec, grad = WiltonInts84.wiltonints(s[1], s[2], s[3], x, Val{1})

    # \int \nabla G_s with G_s = \nabla (1/R + 0.5*γ^2*R) / (4\pi)
    #∫∇G = (-grad[1] - 0.5*γ^2*grad[2]) / (4π)
    ∫∇G = (-grad[1] - 0.5*γ^2*grad[3]) / (4π)

    α = 1 / volume(t) / volume(s) / 4
    for i in 1 : numfunctions(g)
        a = t[i]
        g = (x - a)

        for j in 1 : numfunctions(f)
            b = s[j]

            z[i,j] += ( α * ( (x-b) × g ) ⋅ ∫∇G ) * dx

        end # next j
    end #

end


function select_quadrule()
    try
        Pkg.installed("BogaertInts10")
        info("`BogaertInts10` detected: enhanced quadrature enabled.")
        @eval using BogaertInts10
        @eval include("bogaertints.jl")
        @eval quadrule(op::MaxwellOperator3D, g::RTRefSpace, f::RTRefSpace, i, τ, j, σ, qd) = qrib(op, g, f, i, τ, j, σ, qd)
    catch
        info("Cannot find package `BogaertInts10`. Default quadrature strategy used.")
        @eval quadrule(op::MaxwellOperator3D, g::RTRefSpace, f::RTRefSpace, i, τ, j, σ, qd) = qrdf(op, g, f, i, τ, j, σ, qd)
    end
end

select_quadrule()

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
