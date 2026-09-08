

mutable struct NitscheHH3{T} <: MaxwellOperator3D{T,T}
    gamma::T
end

struct NitscheHH3Reg{T} <: MaxwellOperator3D{T,T}
    gamma::T
end
struct NitscheHH3Sng{T} <: MaxwellOperator3D{T,T}
    gamma::T
end

regularpart(op::NitscheHH3) = NitscheHH3Reg(op.gamma)
singularpart(op::NitscheHH3) = NitscheHH3Sng(op.gamma)

defaultquadstrat(::NitscheHH3, ::LagrangeRefSpace, ::LagrangeRefSpace) = DoubleNumWiltonSauterQStrat(10,8,10,8,3,3,3,3)

function quaddata(operator::NitscheHH3,
    localtestbasis::LagrangeRefSpace,
    localtrialbasis::LagrangeRefSpace,
    testelements, trialelements, qs::DoubleNumWiltonSauterQStrat)

  tqd = quadpoints(localtestbasis,  testelements,  (qs.outer_rule_far,qs.outer_rule_near))
  bqd = quadpoints(x -> localtrialbasis(x), trialelements, (qs.inner_rule_far,qs.inner_rule_near))

  #return QuadData(tqd, bqd)
  return (tpoints=tqd, bpoints=bqd)
end

# Wilton near-field extraction is only implemented for P1 trial functions
_wilton_available(::LagrangeRefSpace, ::LagrangeRefSpace{T,1}) where {T} = true
_wilton_available(::LagrangeRefSpace, ::LagrangeRefSpace) = false

function integrate!(op::NitscheHH3, g::LagrangeRefSpace, f::LagrangeRefSpace, i, τ, j, σ, qd,
        qs::DoubleNumWiltonSauterQStrat,
        out=nothing, test_space=nothing, tptr=nothing, trial_space=nothing, bptr=nothing;
        action::QuadRuleAction=ApplyIntegrate())

    T = eltype(eltype(τ.vertices))
    dmin2 = floatmax(T)
    for t in τ.vertices
        for s in σ.vertices
            dmin2 = min(dmin2, LinearAlgebra.norm_sqr(t-s))
        end
    end

    h2 = volume(σ)
    xtol2 = 0.2 * 0.2
    k2 = abs2(gamma(op))
    if _wilton_available(g, f) && max(dmin2*k2, dmin2/16h2) < xtol2
        qrule = WiltonSERule(
            qd.tpoints[2,i],
            DoubleQuadRule(
                qd.tpoints[2,i],
                qd.bpoints[2,j],),)
        return integrate!(action, out, op, test_space, tptr, τ, trial_space, bptr, σ, qrule)
    end

    qrule = DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j]
    )
    integrate!(action, out, op, test_space, tptr, τ, trial_space, bptr, σ, qrule)
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

function innerintegrals!(op::NitscheHH3Sng, test_neighborhood,
    test_refspace::LagrangeRefSpace,
    trial_refspace::LagrangeRefSpace{T,1},
    test_elements, trial_element, zlocal, quadrature_rule::WiltonSERule, dx) where {T}

    γ = gamma(op)

    s1,s2,s3 = trial_element.vertices

    num_tshapes = numfunctions(test_refspace, domain(test_elements))
    num_bshapes = numfunctions(trial_refspace, domain(trial_element))

    x = cartesian(test_neighborhood)
    tx = normalize(tangents(test_neighborhood, 1))

    # P1 surface curlis constant on the triangle, so one evaluation is sufficient
    q = neighborhood(trial_element, (one(T)/3, one(T)/3))
    fvals = trial_refspace(q)

    scal, _ = WiltonInts84.wiltonints(s1,s2,s3,x,Val{1})
    ∫G = (scal[2] + 0.5*γ^2*scal[4]) / (4π)

    gvals = test_refspace(test_neighborhood)

    for i in 1:num_tshapes
        gx = gvals[i].value
        for j in 1:num_bshapes
            zlocal[i,j] += gx * dot(tx, fvals[j].curl) * ∫G * dx
        end
    end

    return nothing
end

function (igd::Integrand{<:NitscheHH3Reg})(x,y,f,g)
    γ = gamma(igd.operator)

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    green = (expm1(-γ*R) - 0.5*γ^2*R^2) / (4pi*R)

    tx = normalize(tangents(x, 1))

    _integrands(f,g) do fi, gi
        fi.value * dot(tx, green*gi.curl)
    end
end
