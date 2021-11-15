
struct HH3DDoubleLayerNear{K}
  gamma::K
end

function HH3DDoubleLayerNear(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    HH3DDoubleLayerNear(-imag(wavenumber))
  else
    HH3DDoubleLayerNear(wavenumber*im)
  end
end

quaddata(op::HH3DDoubleLayerNear,rs,els) = quadpoints(rs,els,(3,))
quadrule(op::HH3DDoubleLayerNear,refspace,p,y,q,el,qdata) = qdata[1,q]

function kernelvals(op::HH3DDoubleLayerNear,y,p)

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

  ∇G = krn.gradgreen
  nx = krn.nx

  fx = f.value

  ∂G∂n = nx ⋅ ∇G

  return ∂G∂n * fx 
end
