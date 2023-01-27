
struct HH3DDoubleLayerNear{K}
  gamma::K
end

struct HH3DSingleLayerNear{K}
  gamma::K
end

struct HH3DHyperSingularNear{K}
  gamma::K
end

struct HH3DDoubleLayerTransposedNear{K}
  gamma::K
end

function HH3DDoubleLayerNear(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    HH3DDoubleLayerNear(-imag(wavenumber))
  else
    HH3DDoubleLayerNear(wavenumber*im)
  end
end

function HH3DSingleLayerNear(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    HH3DSingleLayerNear(-imag(wavenumber))
  else
    HH3DSingleLayerNear(wavenumber*im)
  end
end

function HH3DHyperSingularNear(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    HH3DHyperSingularNear(-imag(wavenumber))
  else
    HH3DHyperSingularNear(wavenumber*im)
  end
end

function HH3DDoubleLayerTransposedNear(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    HH3DDoubleLayerTransposedNear(-imag(wavenumber))
  else
    HH3DDoubleLayerTranposedNear(wavenumber*im)
  end
end

HH3DNear = Union{HH3DSingleLayerNear, HH3DDoubleLayerNear, HH3DDoubleLayerTransposedNear, HH3DHyperSingularNear}

quaddata(op::HH3DNear,rs,els) = quadpoints(rs,els,(3,))
quadrule(op::HH3DNear,refspace,p,y,q,el,qdata) = qdata[1,q]

function kernelvals(op::HH3DNear,y,p)

  γ = op.gamma
  x = cartesian(p)
  r = y - x
  R = norm(r)
  γR = γ*R

  inv_R = 1/R

  expn = exp(-γR)
  green = expn * inv_R * inv_4pi
  gradgreen = -(γ + inv_R) * green * inv_R * r

  nx = normal(p)

  (;γ, r, R, green, gradgreen, nx)
end

function integrand(op::HH3DDoubleLayerNear,krn,y,f,p)

  ∇G = -krn.gradgreen
  nx = krn.nx

  fx = f.value

  ∂G∂n = nx ⋅ ∇G

  return ∂G∂n * fx 
end

function integrand(op::HH3DSingleLayerNear, krn, y, f, p)
  G = krn.green
  fx = f.value

  return G * fx
end

function integrand(op::HH3DDoubleLayerTransposedNear, krn, y, f, p)
  ∇G = krn.gradgreen
  fx = f.value

  return ∇G * fx
end

function integrand(op::HH3DHyperSingularNear, krn, y, f, p)
  G = krn.green
  nx = krn.nx
  γ = krn.γ
  invR = 1/krn.R
  fx = f.value
  r = krn.r

  # returns ∇ᵣn̂'∇ᵣ'G - strong singularity, probably makes problems for small distances R
  return (nx * G * invR * (γ + invR) - r * dot(nx,r) * G * invR^2 * (γ^2 + 3 * γ * invR + 3 * invR^2)) * fx
end
