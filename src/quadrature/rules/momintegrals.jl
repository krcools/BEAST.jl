function momintegrals!(out, op,
    test_functions, test_cellptr, test_chart,
    trial_functions, trial_cellptr, trial_chart,
    quadrule)

    local_test_space = refspace(test_functions)
    local_trial_space = refspace(trial_functions)

    momintegrals!(op,
        local_test_space, local_trial_space,
        test_chart, trial_chart,
        out, quadrule)
end