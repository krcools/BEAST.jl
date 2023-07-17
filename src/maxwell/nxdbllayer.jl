

mutable struct DoubleLayerRotatedMW3D{T} <: MaxwellOperator3D
    "im times the wavenumber"
    gamma::T
end

LinearAlgebra.cross(::NormalVector, a::MWDoubleLayer3D) = DoubleLayerRotatedMW3D(a.gamma)

# defaultquadstrat(::DoubleLayerRotatedMW3D, tfs, bfs) = DoubleNumQStrat(2,3)

function quaddata(operator::DoubleLayerRotatedMW3D,
        local_test_basis::LinearRefSpaceTriangle, local_trial_basis::LinearRefSpaceTriangle,
        test_elements, trial_elements, qs::DoubleNumQStrat)

    test_quad_data  = quadpoints(local_test_basis,  test_elements,  (qs.outer_rule,))
    trial_quad_data = quadpoints(local_trial_basis, trial_elements, (qs.inner_rule,))

    return test_quad_data, trial_quad_data
end


function quadrule(operator::DoubleLayerRotatedMW3D,
        local_test_basis::LinearRefSpaceTriangle, local_trial_basis::LinearRefSpaceTriangle,
        test_id, test_element, trial_id, trial_element,
        quad_data, qs::DoubleNumQStrat)

    test_quad_rules  = quad_data[1]
    trial_quad_rules = quad_data[2]

    DoubleQuadRule(
        test_quad_rules[1,test_id],
        trial_quad_rules[1,trial_id]
    )
end

# TODO: this method needs to go and dispatch needs to be dealt with using the defaultquadstrat mechanism
function quadrule(op::DoubleLayerRotatedMW3D, g::RTRefSpace, f::RTRefSpace,  i, τ, j, σ, qd,
    qs::DoubleNumWiltonSauterQStrat)
  
    hits = 0
    dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))
    dmin2 = floatmax(eltype(eltype(τ.vertices)))
    for t in τ.vertices
        for s in σ.vertices
            d2 = LinearAlgebra.norm_sqr(t-s)
            dmin2 = min(dmin2, d2)
            hits += (d2 < dtol)
        end
    end
  
    hits == 3 && return SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
    hits == 2 && return SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
    hits == 1 && return SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])
  
    h2 = volume(σ)
    xtol2 = 0.2 * 0.2
    k2 = abs2(op.gamma)
    return DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j],)
end

# function integrand(op::DoubleLayerRotatedMW3D, kernel_vals, test_vals, test_nbd, trial_vals, trial_nbd)

#     n = normal(test_nbd)
#     g = test_vals[1]
#     f = trial_vals[1]
#     ∇G = kernel_vals.gradgreen

#     return g ⋅ (n × (∇G × f))
# end

# function (igd::Integrand{<:DoubleLayerRotatedMW3D})(u,v)

#     x = neighborhood(igd.test_chart,u)
#     y = neighborhood(igd.trial_chart,v)
#     j = jacobian(x) * jacobian(y)
#     nx = normal(x)

#     r = cartesian(x) - cartesian(y)
#     R = norm(r)
#     iR = 1/R
#     γ = igd.operator.gamma
#     G = exp(-γ*R)/(4π*R)
#     K = -(γ + iR) * G * (iR * r)

#     f = igd.local_test_space(x)
#     g = igd.local_trial_space(y)

#     fvalue = getvalue(f)
#     gvalue = getvalue(g)

#     jKg = cross.(Ref(K), j*gvalue)
#     jnxKg = cross.(Ref(nx), jKg)
#     return _krondot(fvalue, jnxKg)
# end


function (igd::Integrand{<:DoubleLayerRotatedMW3D})(x,y,f,g)

    # x = neighborhood(igd.test_chart,u)
    # y = neighborhood(igd.trial_chart,v)
   # j = jacobian(x) * jacobian(y)
    nx = normal(x)

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    γ = igd.operator.gamma
    G = exp(-γ*R)/(4π*R)
    K = -(γ + iR) * G * (iR * r)

    # f = igd.local_test_space(x)
    # g = igd.local_trial_space(y)

    fvalue = getvalue(f)
    gvalue = getvalue(g)

    #jKg = cross.(Ref(K), j*gvalue)
    jKg = cross.(Ref(K), gvalue)
    jnxKg = cross.(Ref(nx), jKg)
    return _krondot(fvalue, jnxKg)
end
