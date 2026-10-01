struct TestRefinesTrialQRule{S}
    conforming_qstrat::S
end

quadraturebuffer(qr::TestRefinesTrialQRule, test_space, trial_space) =
    quadraturebuffer(qr.conforming_qstrat, test_space, trial_space)

function integrate!(out, op,
    test_functions::Space, test_cell, test_chart,
    trial_functions::Space, trial_cell, trial_chart,
    qr::TestRefinesTrialQRule)

    return integrate!(
        out, op, test_functions, test_cell, test_chart,
        trial_functions, trial_cell, trial_chart,
        qr, quadraturebuffer(qr, test_functions, trial_functions))
end

function integrate!(out, op,
    test_functions::Space, test_cell, test_chart,
    trial_functions::Space, trial_cell, trial_chart,
    qr::TestRefinesTrialQRule, qbuffer)

    test_local_space = refspace(test_functions)
    trial_local_space = refspace(trial_functions)

    test_mesh = geometry(test_functions)
    trial_mesh = geometry(trial_functions)

    tdom = domain(test_chart)
    bdom = domain(trial_chart)

    num_tshapes = numfunctions(test_local_space, tdom)
    num_bshapes = numfunctions(trial_local_space, bdom)

    parent_mesh = CompScienceMeshes.parent(test_mesh)
    trial_charts = [chart(test_mesh, p) for p in CompScienceMeshes.children(parent_mesh, trial_cell)]

    trial_overlaps = map(trial_charts) do chart
        simplex(map(v -> carttobary(trial_chart, v), chart.vertices))
    end

    quadstrat = qr.conforming_qstrat
    qd = quaddata(op, test_local_space, trial_local_space,
        [test_chart], trial_charts, quadstrat)

    zlocal = zero(out)
    Q = zeros(coordtype(trial_chart), num_bshapes, num_bshapes)
    qaction = ApplyIntegrate(qbuffer)
    for (q,chart) in enumerate(trial_charts)
        restrict!(Q, trial_local_space, trial_chart, chart, trial_overlaps[q])

        fill!(zlocal, 0)
        integrate!(op, test_local_space, trial_local_space,
            1, test_chart, q, chart, qd, quadstrat,
            zlocal, test_functions, nothing, trial_functions, nothing;
            action=qaction)

        for j in 1:num_bshapes
            for i in 1:num_tshapes
                for k in 1:size(Q, 2)
                    out[i,j] += zlocal[i,k] * Q[j,k]
end end end end end
