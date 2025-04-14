struct TestInBaryRefOfTrialQRule{S}
    conforming_qstrat::S
end

function BEAST.momintegrals!(out, op,
    test_functions, test_cell, test_chart,
    trial_functions, trial_cell, trial_chart,
    qr::TestInBaryRefOfTrialQRule)

    test_local_space = refspace(test_functions)
    trial_local_space = refspace(trial_functions)

    num_tshapes = numfunctions(test_local_space, domain(test_chart))
    num_bshapes = numfunctions(trial_local_space, domain(trial_chart))

    T = coordtype(test_chart)
    z, u, h, t = zero(T), one(T), T(1//2), T(1//3)

    c = CompScienceMeshes.point(T, t, t)
    v = (
        CompScienceMeshes.point(T, u, z),
        CompScienceMeshes.point(T, z, u),
        CompScienceMeshes.point(T, z, z))
    e = (
        CompScienceMeshes.point(T, z, h),
        CompScienceMeshes.point(T, h, z),
        CompScienceMeshes.point(T, h, h))

    X = (
        CompScienceMeshes.simplex(v[1], e[3], c),
        CompScienceMeshes.simplex(v[2], c, e[3]),
        CompScienceMeshes.simplex(v[2], e[1], c),
        CompScienceMeshes.simplex(v[3], c, e[1]),
        CompScienceMeshes.simplex(v[3], e[2], c),
        CompScienceMeshes.simplex(v[1], c, e[2]))

    C = CompScienceMeshes.cartesian(trial_chart, c)
    V = (
        CompScienceMeshes.cartesian(trial_chart, v[1]),
        CompScienceMeshes.cartesian(trial_chart, v[2]),
        CompScienceMeshes.cartesian(trial_chart, v[3]))
    E = (
        CompScienceMeshes.cartesian(trial_chart, e[1]),
        CompScienceMeshes.cartesian(trial_chart, e[2]),
        CompScienceMeshes.cartesian(trial_chart, e[3]))

    trial_charts = (
        CompScienceMeshes.simplex(V[1], E[3], C),
        CompScienceMeshes.simplex(V[2], C, E[3]),
        CompScienceMeshes.simplex(V[2], E[1], C),
        CompScienceMeshes.simplex(V[3], C, E[1]),
        CompScienceMeshes.simplex(V[3], E[2], C),
        CompScienceMeshes.simplex(V[1], C, E[2]))

    quadstrat = qr.conforming_qstrat
    qd = BEAST.quaddata(op, test_local_space, trial_local_space,
        (test_chart,), trial_charts, quadstrat)

    Q = zeros(T, num_tshapes, num_tshapes)
    out1 = zero(out)
    for (q,chart) in enumerate(trial_charts)
        qr1 = BEAST.quadrule(op, test_local_space, trial_local_space,
            1, test_chart, q ,chart, qd, quadstrat)
            
        BEAST.restrict!(Q, trial_local_space, trial_chart, chart, X[q])

        fill!(out1, 0)
        BEAST.momintegrals!(out1, op,
            test_functions, nothing, test_chart,
            trial_functions, nothing, chart, qr1)

        for j in 1:num_bshapes
            for i in 1:num_tshapes
                for k in 1:size(Q, 2)
                    out[i,j] += out1[i,k] * Q[j,k]
end end end end end
