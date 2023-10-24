struct SymmetricQuadStrat{S}
    quadstrat::S
end

struct SymmetricQuadRule{R1,R2}
    quadrule1::R1
    quadrule2::R2
end

function quaddata(op,
    test_local_space, trial_local_space, time_local_space,
    test_charts, trial_charts, time_charts,
    quadstrat::SymmetricQuadStrat)

    qd1 = quaddata(op,
        test_local_space, trial_local_space, time_local_space,
        test_charts, trial_charts, time_charts,
        quadstrat.quadstrat)
    qd2 = quaddata(op,
        trial_local_space, test_local_space, time_local_space,
        trial_charts, test_charts, time_charts,
        quadstrat.quadstrat)

    return qd1, qd2
end

function quadrule(op,
    test_local_space, trial_local_space, time_local_space,
    p, test_chart, q, trial_chart, r, time_chart,
    qd, quadstrat::SymmetricQuadStrat)

    qd1 = qd[1]
    qd2 = qd[2]

    qr1 = quadrule(op, test_local_space, trial_local_space, time_local_space,
        p, test_chart, q, trial_chart, r, time_chart,
        qd1, quadstrat.quadstrat)
    qr2 = quadrule(op, trial_local_space, test_local_space, time_local_space,
        q, trial_chart, p, test_chart, r, time_chart,
        qd2, quadstrat.quadstrat)

    return SymmetricQuadRule(qr1, qr2)
end

function momintegrals!(z, op, U, V, W, τ, σ, ι, qr::SymmetricQuadRule)

    qr1 = qr.quadrule1
    qr2 = qr.quadrule2

    z1 = zero(z)
    z2 = zero(z)

    momintegrals!(z1, op, U, V, W, τ, σ, ι, qr1)
    momintegrals!(z2, op, V, U, W, σ, τ, ι, qr2)
    z2 = permutedims(z2, (2,1,3))
    z .+= (z1+z2)/2
    return nothing
end