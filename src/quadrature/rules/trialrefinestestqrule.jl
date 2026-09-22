struct TrialRefinesTestQRule{S}
    conforming_qstrat::S
end

quadraturebuffer(qr::TrialRefinesTestQRule, test_space, trial_space) =
    quadraturebuffer(qr.conforming_qstrat, test_space, trial_space)

function integrate!(out, op,
    test_functions::Space, test_cell, test_chart,
    trial_functions::Space, trial_cell, trial_chart,
    qr::TrialRefinesTestQRule)

    return integrate!(
        out, op, test_functions, test_cell, test_chart,
        trial_functions, trial_cell, trial_chart,
        qr, quadraturebuffer(qr, test_functions, trial_functions))
end

function integrate!(out, op,
    test_functions::Space, test_cell, test_chart,
    trial_functions::Space, trial_cell, trial_chart,
    qr::TrialRefinesTestQRule, qbuffer)

    test_local_space = refspace(test_functions)
    trial_local_space = refspace(trial_functions)

    test_mesh = geometry(test_functions)
    trial_mesh = geometry(trial_functions)

    tdom = domain(chart(test_mesh, first(test_mesh)))
    bdom = domain(chart(trial_mesh, first(trial_mesh)))

    num_tshapes = numfunctions(test_local_space, tdom)
    num_bshapes = numfunctions(trial_local_space, bdom)

    parent_mesh = CompScienceMeshes.parent(trial_mesh)
    test_charts = [chart(trial_mesh, p) for p in CompScienceMeshes.children(parent_mesh, test_cell)]

    quadstrat = qr.conforming_qstrat
    qd = quaddata(op, test_local_space, trial_local_space,
        test_charts, [trial_chart], quadstrat)
    qaction = ApplyIntegrate(qbuffer)

    for (p,chart) in enumerate(test_charts)
        Q = restrict(test_local_space, test_chart, chart)
        zlocal = zero(out)
        integrate!(op, test_local_space, trial_local_space,
            p, chart, 1, trial_chart, qd, quadstrat,
            zlocal, test_functions, nothing, trial_functions, trial_cell;
            action=qaction)

        for j in 1:num_bshapes
            for i in 1:num_tshapes
                for k in 1:size(Q, 2)
                    out[i,j] += Q[i,k] * zlocal[k,j]
end end end end end
