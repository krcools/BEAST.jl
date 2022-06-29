

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

function integrand(op::DoubleLayerRotatedMW3D, kernel_vals, test_vals, test_nbd, trial_vals, trial_nbd)

    n = normal(test_nbd)
    g = test_vals[1]
    f = trial_vals[1]
    ∇G = kernel_vals.gradgreen

    return g ⋅ (n × (∇G × f))
end
