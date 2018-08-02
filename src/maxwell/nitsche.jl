

"""
Describe a single layer operator from the surface to a line.


```math
<v, Su> = ∫_γ dx v(x) ∫_Γ dy \frac{e^{-ikR}}{4πR} u(y)
```
"""
mutable struct SingleLayerTrace{T} <: MaxwellOperator3D
    gamma::T
end

function quaddata(operator::SingleLayerTrace,
  localtestbasis::LagrangeRefSpace,
  localtrialbasis::LagrangeRefSpace,
  testelements,
  trialelements)

  tqd = quadpoints(localtestbasis,  testelements,  (10,))
  bqd = quadpoints(localtrialbasis, trialelements, (8,))

  #return QuadData(tqd, bqd)
  return (tpoints=tqd, bpoints=bqd)
end

# Use numerical quadrature for now
# Note: basis integral is over triangle, test over line
function quadrule(op::SingleLayerTrace, g::LagrangeRefSpace, f::LagrangeRefSpace, i, τ, j, σ, qd)
    DoubleQuadStrategy(
        qd.tpoints[1,i],
        qd.bpoints[1,j]
    )
end

integrand(op::SingleLayerTrace, kernel, g, τ, f, σ) = f[1]*g[1]*kernel.green
