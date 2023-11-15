

mutable struct NitscheHH3{T} <: MaxwellOperator3D{T,T}
    gamma::T
end

defaultquadstrat(::NitscheHH3, ::LagrangeRefSpace, ::LagrangeRefSpace) = DoubleNumWiltonSauterQStrat(10,8,10,8,3,3,3,3)

function quaddata(operator::NitscheHH3,
    localtestbasis::LagrangeRefSpace,
    localtrialbasis::LagrangeRefSpace,
    testelements, trialelements, qs::DoubleNumWiltonSauterQStrat)

  tqd = quadpoints(localtestbasis,  testelements,  (qs.outer_rule_far,))
  bqd = quadpoints(x -> localtrialbasis(x), trialelements, (qs.inner_rule_far,))

  #return QuadData(tqd, bqd)
  return (tpoints=tqd, bpoints=bqd)
end

function quadrule(op::NitscheHH3, g::LagrangeRefSpace, f::LagrangeRefSpace, i, τ, j, σ, qd,
        qs::DoubleNumWiltonSauterQStrat)
        
    DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j]
    )
end


struct KernelValsMaxwell3D{T,U,P,Q}
    "gamma = im * wavenumber"
    gamma::U
    vect::P
    dist::T
    green::U
    gradgreen::Q
end

const inv_4pi = 1/(4pi)
function kernelvals(biop::MaxwellOperator3D, p, q)

    γ = gamma(biop)
    r = cartesian(p) - cartesian(q)
    T = eltype(r)
    R = norm(r)
    γR = γ*R

    inv_R = 1/R

    expn = exp(-γR)
    green = expn * inv_R * T(inv_4pi)
    gradgreen = -(γ + inv_R) * green * inv_R * r

    KernelValsMaxwell3D(γ, r, R, green, gradgreen)
end

function integrand(op::NitscheHH3, kernel, test_vals, test_point, trial_vals, trial_point)
    Gxy = kernel.green
    @assert length(test_point.patch.tangents) == 1
    tx = normalize(tangents(test_point, 1))
    gx = test_vals[1]
    curlfy = trial_vals[2]
    return gx*dot(tx, Gxy * curlfy)
end
