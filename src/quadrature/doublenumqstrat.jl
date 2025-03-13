struct DoubleNumQStrat{R} <: AbstractQuadStrat
    outer_rule::R
    inner_rule::R
end

function quaddata(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_elements, trial_elements, qs::DoubleNumQStrat)

    # local_test_basis = refspace(test_basis)
    # local_trial_basis = refspace(trial_basis)

    test_quad_data  = quadpoints(local_test_basis,  test_elements,  (qs.outer_rule,))
    trial_quad_data = quadpoints(local_trial_basis, trial_elements, (qs.inner_rule,))

    return test_quad_data, trial_quad_data
end


function quadrule(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,
    quad_data, qs::DoubleNumQStrat)

    test_quad_rules  = quad_data[1]
    trial_quad_rules = quad_data[2]

    DoubleQuadRule(
        test_quad_rules[1,test_id],
        trial_quad_rules[1,trial_id]
    )
end