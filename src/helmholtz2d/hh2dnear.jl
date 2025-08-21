
struct HH2DDoubleLayerNear{K}
  gamma::K
end

struct HH2DSingleLayerNear{K}
  gamma::K
end

struct HH2DHyperSingularNear{K}
  gamma::K
end

struct HH2DDoubleLayerTransposedNear{K}
  gamma::K
end

function HH2DDoubleLayerNear(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    HH2DDoubleLayerNear(-imag(wavenumber))
  else
    HH2DDoubleLayerNear(wavenumber*im)
  end
end

function HH2DSingleLayerNear(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    HH2DSingleLayerNear(-imag(wavenumber))
  else
    HH2DSingleLayerNear(wavenumber*im)
  end
end

function HH2DHyperSingularNear(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    HH2DHyperSingularNear(-imag(wavenumber))
  else
    HH2DHyperSingularNear(wavenumber*im)
  end
end

function HH2DDoubleLayerTransposedNear(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    HH2DDoubleLayerTransposedNear(-imag(wavenumber))
  else
    HH2DDoubleLayerTranposedNear(wavenumber*im)
  end
end

HH2DNear = Union{HH2DSingleLayerNear, HH2DDoubleLayerNear, HH2DDoubleLayerTransposedNear, HH2DHyperSingularNear}
defaultquadstrat(op::HH2DNear, basis) = SingleNumQStrat(20)
quaddata(op::HH2DNear,rs,els,qs::SingleNumQStrat) = quadpoints(rs,els,(qs.quad_rule,))
quadrule(op::HH2DNear,refspace,p,y,q,el,qdata,qs::SingleNumQStrat) = qdata[1,q]

# TODO: This function is useful for changing the gradient to a curl
# However, it is a case of type piracy.
# (Evidently, the case is not supported by cross function due to mismatch of length)

@inline function LinearAlgebra.cross(z::SVector{3,T1}, v::SVector{2,T2}) where {T1,T2}
    z == ẑ || throw(ArgumentError("Only ẑ = (0,0,1) supported"))
    return SVector{2,T2}(-v[2],  v[1])
end

@inline function LinearAlgebra.cross(v::SVector{2,T2}, z::SVector{3,T1}) where {T1,T2}
    z == ẑ || throw(ArgumentError("Only ẑ = (0,0,1) supported"))
    return SVector{2,T2}( v[2], -v[1])
end

# WARNING: Matching the HH3D we use x for the sources
# TODO: Discuss if we change this?
function kernelvals(op::HH2DNear,y,p)

    γ = op.gamma
    if iszero(real(γ))
        k = imag(γ)
    else
        k = -im*γ
    end
    x = cartesian(p)
    r = y - x
    R = norm(r)

    kR = k * R
    hankels = hankelh2.([0 1], kR)
    #green = - op.alpha * im / 4 * hankels[1]
    #gradgreen = op.alpha * k * im / 4 * hankels[2] * r / R
    green = - im / 4 * hankels[1]
    gradgreen =  k * im / 4 * hankels[2] * r / R

    ddhankel = hankels[1] - 1/(kR) * hankels[2]

    #txty = dot(normal(tgeo), normal(bgeo))

    nx = normal(p)

    # Needed for hypersingular operator
    gradnxgradgreen = ddhankel * (k^2/R^2 * dot(r, nx) * r)
    gradnxgradgreen += k * hankels[2] * (nx / R - r/R^3 * dot(r, nx)) 
    gradnxgradgreen *= - im / 4

    (;γ, r, R, green, gradgreen, nx, gradnxgradgreen)
end

function integrand(op::HH2DDoubleLayerNear,krn,y,f,p)

  ∇G = -krn.gradgreen
  nx = krn.nx

  fx = f.value

  ∂G∂n = nx ⋅ ∇G

  return ∂G∂n * fx 
end

function integrand(op::HH2DSingleLayerNear, krn, y, f, p)
  G = krn.green
  fx = f.value

  return G * fx
end

function integrand(op::HH2DDoubleLayerTransposedNear, krn, y, f, p)
  ∇G = krn.gradgreen
  fx = f.value

  return ∇G * fx
end

function integrand(op::HH2DHyperSingularNear, krn, y, f, p)
  ∇nₓ∇ₓG = krn.gradnxgradgreen
  fx = f.value

  return ∇nₓ∇ₓG * fx
end
