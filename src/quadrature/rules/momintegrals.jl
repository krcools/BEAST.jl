function integrate!(out, op,
    test_functions, test_cellptr, test_chart,
    trial_functions, trial_cellptr, trial_chart,
    quadrule)

    return integrate!(
        out, op, test_functions, test_cellptr, test_chart,
        trial_functions, trial_cellptr, trial_chart,
        quadrule, quadraturebuffer(quadrule, test_functions, trial_functions))
end

function integrate!(out, op,
    test_functions, test_cellptr, test_chart,
    trial_functions, trial_cellptr, trial_chart,
    quadrule, qbuffer)

    local_test_space = refspace(test_functions)
    local_trial_space = refspace(trial_functions)

    integrate!(op,
        local_test_space, local_trial_space,
        test_chart, trial_chart,
        out, quadrule, qbuffer)
end
