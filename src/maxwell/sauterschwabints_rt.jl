struct MWSL3DIntegrand{C,O,L}
    test_triangular_element::C
    trial_triangular_element::C
    op::O
    test_local_space::L
    trial_local_space::L
end

function (igd::MWSL3DIntegrand)(u,v)
        α = igd.op.α
        β = igd.op.β
        γ = igd.op.gamma

        x = neighborhood(igd.test_triangular_element,u)
        y = neighborhood(igd.trial_triangular_element,v)

        r = cartesian(x) - cartesian(y)
        R = norm(r)
        G = exp(-γ*R)/(4π*R)

        f = igd.test_local_space(x)
        g = igd.trial_local_space(y)

        j = jacobian(x) * jacobian(y)

        αjG = α*j*G
        βjG = β*j*G

        G = @SVector [αjG*g[1].value, αjG*g[2].value, αjG*g[3].value]
        H = @SVector [βjG*g[1].divergence, βjG*g[2].divergence, βjG*g[3].divergence]

        SMatrix{3,3}((
            dot(f[1].value,G[1])+f[1].divergence*H[1],
            dot(f[2].value,G[1])+f[2].divergence*H[1],
            dot(f[3].value,G[1])+f[3].divergence*H[1],
            dot(f[1].value,G[2])+f[1].divergence*H[2],
            dot(f[2].value,G[2])+f[2].divergence*H[2],
            dot(f[3].value,G[2])+f[3].divergence*H[2],
            dot(f[1].value,G[3])+f[2].divergence*H[3],
            dot(f[2].value,G[3])+f[2].divergence*H[3],
            dot(f[3].value,G[3])+f[3].divergence*H[3]))
end

struct MWDL3DIntegrand{C,O,L}
        test_triangular_element::C
        trial_triangular_element::C
        op::O
        test_local_space::L
        trial_local_space::L
end

function (igd::MWDL3DIntegrand)(u,v)

        γ = igd.op.gamma

        x = neighborhood(igd.test_triangular_element,u)
        y = neighborhood(igd.trial_triangular_element,v)

        r = cartesian(x) - cartesian(y)
        R = norm(r)
        G = exp(-γ*R)/(4π*R)

        f = igd.test_local_space(x)
        g = igd.trial_local_space(y)

        j = jacobian(x) * jacobian(y)

        GG = -(γ + 1/R) * G / R * r
        T = @SMatrix [
        0 -GG[3] GG[2]
        GG[3] 0 -GG[1]
        -GG[2] GG[1] 0 ]

        G = @SVector [j*T*g[1].value, j*T*g[2].value, j*T*g[3].value]
        SMatrix{3,3}((
            dot(f[1].value, G[1]),
            dot(f[2].value, G[1]),
            dot(f[3].value, G[1]),
            dot(f[1].value, G[2]),
            dot(f[2].value, G[2]),
            dot(f[3].value, G[2]),
            dot(f[1].value, G[3]),
            dot(f[2].value, G[3]),
            dot(f[3].value, G[3])))
end

kernel_in_bary(op::MWSingleLayer3D,
    test_local_space::RTRefSpace, trial_local_space::RTRefSpace,
    test_chart, trial_chart) = MWSL3DIntegrand(
        test_chart, trial_chart, op, test_local_space, trial_local_space)

kernel_in_bary(op::MWDoubleLayer3D,
    test_local_space::RTRefSpace, trial_local_space::RTRefSpace,
    test_chart, trial_chart) = MWDL3DIntegrand(
        test_chart, trial_chart, op, test_local_space, trial_local_space)

const MWOperator3D = Union{MWSingleLayer3D, MWDoubleLayer3D}
function momintegrals!(op::MWOperator3D,
    test_local_space::RTRefSpace, trial_local_space::RTRefSpace,
    test_triangular_element, trial_triangular_element, out, strat::SauterSchwabStrategy)

    I, J, K, L = SauterSchwabQuadrature.reorder(
        test_triangular_element.vertices,
        trial_triangular_element.vertices, strat)

    test_triangular_element  = simplex(
        test_triangular_element.vertices[I[1]],
        test_triangular_element.vertices[I[2]],
        test_triangular_element.vertices[I[3]])

    trial_triangular_element = simplex(
        trial_triangular_element.vertices[J[1]],
        trial_triangular_element.vertices[J[2]],
        trial_triangular_element.vertices[J[3]])

    # igd = MWSL3DIntegrand(test_triangular_element, trial_triangular_element,
    #     op, test_local_space, trial_local_space)
    igd = kernel_in_bary(op, test_local_space, trial_local_space,
        test_triangular_element, trial_triangular_element)
    G = sauterschwab_parameterized(igd, strat)
    for j ∈ 1:3, i ∈ 1:3
        out[i,j] += G[K[i],L[j]]
    end

    nothing
end

function momintegrals_nested!(op::MWOperator3D,
    test_local_space::RTRefSpace, trial_local_space::RTRefSpace,
    test_chart, trial_chart, out, strat::SauterSchwabStrategy)

    # 1. Refine the trial_chart
    p1, p2, p3 = trial_chart.vertices
    e1 = cartesian(neighborhood(trial_chart, (0,1/2)))

    e2 = cartesian(neighborhood(trial_chart, (1/2,0)))
    e3 = cartesian(neighborhood(trial_chart, (1/2,1/2)))
    ct = cartesian(center(trial_chart))
    refined_trial_chart = [
        simplex(ct, p1, e3),
        simplex(ct, e3, p2),
        simplex(ct, p2, e1),
        simplex(ct, e1, p3),
        simplex(ct, p3, e2),
        simplex(ct, e2, p1)]

    qd = quaddata(op, test_local_space, trial_local_space,
        [test_chart], refined_trial_chart)

    for (q,chart) in enumerate(refined_trial_chart)
        qr = quadrule(op, test_local_space, trial_local_space,
            1, test_chart, q ,chart, qd)

        Q = restrict(trial_local_space, trial_chart, chart)
        zlocal = zero(out)
        momintegrals!(op, test_local_space, trial_local_space,
            test_chart, chart, zlocal, qr)

        for j in 1:3
            for i in 1:3
                for k in 1:3
                out[i,j] += zlocal[i,k] * Q[j,k]
        end end end
    end
end

function momintegrals_nested!(op::IntegralOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_chart, trial_chart, out, strat)

    momintegrals!(op, test_local_space, trial_local_space,
        test_chart, trial_chart, out, strat)
end
