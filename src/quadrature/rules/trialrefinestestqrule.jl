struct TrialRefinesTestQRule{S}
    conforming_qstrat::S
end

function momintegrals!(out, op,
    test_functions::Space, test_cell, test_chart,
    trial_functions::Space, trial_cell, trial_chart,
    qr::TrialRefinesTestQRule)

    test_local_space = refspace(test_functions)
    trial_local_space = refspace(trial_functions)

    test_mesh = geometry(test_functions)
    trial_mesh = geometry(trial_functions)

    parent_mesh = CompScienceMeshes.parent(trial_mesh)
    test_charts = [chart(trial_mesh, p) for p in CompScienceMeshes.children(parent_mesh, test_cell)]

    quadstrat = qr.conforming_qstrat
    qd = quaddata(op, test_local_space, trial_local_space,
        test_charts, [trial_chart], quadstrat)

    for (p,chart) in enumerate(test_charts)
        qr = quadrule(op, test_local_space, trial_local_space,
            p, chart, 1, trial_chart, qd, quadstrat)

        Q = restrict(test_local_space, test_chart, chart)
        zlocal = zero(out)
        momintegrals!(zlocal, op,
            test_functions, nothing, chart,
            trial_functions, trial_cell, trial_chart, qr)

        for j in 1:numfunctions(trial_local_space)
            for i in 1:numfunctions(test_local_space)
                for k in 1:size(Q, 2)
                    out[i,j] += Q[i,k] * zlocal[k,j]
end end end end end