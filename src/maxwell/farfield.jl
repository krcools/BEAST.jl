

"""
Operator to compute the far field of a current distribution. In particular, given the current distribution ``j`` this operator allows for the computation of

```math
A j = n × ∫_Γ e^{γ ̂x ⋅ y} dy
```

where ``̂x`` is the unit vector in the direction of observation. Note that the assembly routing expects the observation directions to be normalised by the caller.
"""
struct MWFarField3D{K}
  gamma::K
end

function MWFarField3D(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    MWFarField3D(-imag(wavenumber))
  else
    MWFarField3D(wavenumber*im)
  end
end

quaddata(op::MWFarField3D,rs,els) = quadpoints(rs,els,(3,))
quadrule(op::MWFarField3D,refspace,p,y,q,el,qdata) = qdata[1,q]

kernelvals(op::MWFarField3D,y,p) = exp(op.gamma*dot(y,cartesian(p)))
function integrand(op::MWFarField3D,krn,y,f,p)
    (y × (krn * f[1])) × y
end
