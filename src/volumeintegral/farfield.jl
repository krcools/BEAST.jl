"""
Operator to compute the far field of a current distribution. In particular, given the current distribution ``j`` this operator allows for the computation of

```math
A j = n × ∫_Ω j e^{γ x ⋅ y} dy
```

where ``x`` is the unit vector in the direction of observation. Note that the assembly routing expects the observation directions to be normalised by the caller.
"""
struct VIEFarField3D{K,P}
  gamma::K
  tau::P
end

function VIEFarField3D(;wavenumber=error("wavenumber is a required argument"), tau=nothing)
  gamma = nothing

  if iszero(real(wavenumber))
    gamma = -imag(wavenumber)
  else
    gamma = wavenumber*im
  end

  tau == nothing && (tau = x->1.0)

  VIEFarField3D(gamma, tau)
end

quaddata(op::VIEFarField3D,rs,els) = quadpoints(rs,els,(3,))
quadrule(op::VIEFarField3D,refspace,p,y,q,el,qdata) = qdata[1,q]

kernelvals(op::VIEFarField3D,y,p) = exp(op.gamma*dot(y,cartesian(p)))
function integrand(op::VIEFarField3D,krn,y,f,p)
    (y × (op.tau(cartesian(p)) * krn * f[1])) × y
end
