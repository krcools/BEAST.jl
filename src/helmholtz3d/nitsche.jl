

mutable struct NitscheHH3{T} <: MaxwellOperator3D
    gamma::T
end

function quaddata(operator::NitscheHH3,
    localtestbasis::LagrangeRefSpace,
    localtrialbasis::LagrangeRefSpace,
    testelements, trialelements)

  tqd = quadpoints(localtestbasis,  testelements,  (10,))
  bqd = quadpoints(x -> localtrialbasis(x, Val{:withcurl}), trialelements, (8,))

  #return QuadData(tqd, bqd)
  return (tpoints=tqd, bpoints=bqd)
end

function quadrule(op::NitscheHH3, g::LagrangeRefSpace, f::LagrangeRefSpace, i, τ, j, σ, qd)
    DoubleQuadStrategy(
        qd.tpoints[1,i],
        qd.bpoints[1,j]
    )
end

function integrand(op::NitscheHH3, kernel, test_vals, test_point, trial_vals, trial_point)
    Gxy = kernel.green
    @assert length(test_point.patch.tangents) == 1
    tx = normalize(tangents(test_point, 1))
    gx = test_vals[1]
    curlfy = trial_vals[2]
    return gx*dot(tx, Gxy * curlfy)
end
