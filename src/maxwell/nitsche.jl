

"""
Describe a single layer operator from the surface to a line.


```math
<v, Su> = ∫_γ dx v(x) ∫_Γ dy \frac{e^{-ikR}}{4πR} u(y)
```
"""
mutable struct SingleLayerTrace{T} <: MaxwellOperator3D{T,T}
    gamma::T
end

defaultquadstrat(::SingleLayerTrace, ::LagrangeRefSpace, ::LagrangeRefSpace) = DoubleNumWiltonSauterQStrat(10,8,10,8,3,3,3,3)
function quaddata(operator::SingleLayerTrace,
  localtestbasis::LagrangeRefSpace,
  localtrialbasis::LagrangeRefSpace,
  testelements, trialelements, qs::DoubleNumWiltonSauterQStrat)

  tqd = quadpoints(localtestbasis,  testelements,  (qs.outer_rule_far,))
  bqd = quadpoints(localtrialbasis, trialelements, (qs.outer_rule_near,))

  #return QuadData(tqd, bqd)
  return (tpoints=tqd, bpoints=bqd)
end

# Use numerical quadrature for now
# Note: basis integral is over triangle, test over line
function quadrule(op::SingleLayerTrace, g::LagrangeRefSpace, f::LagrangeRefSpace, i, τ, j, σ, qd,
        qs::DoubleNumWiltonSauterQStrat)
        
    DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j]
    )
end

integrand(op::SingleLayerTrace, kernel, g, τ, f, σ) = f[1]*g[1]*kernel.green
