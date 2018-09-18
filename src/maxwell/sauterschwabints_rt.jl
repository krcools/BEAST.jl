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

        G = @SVector [αjG*g[i].value for i in 1:3]
        H = @SVector [βjG*g[i].divergence for i in 1:3]

        @SMatrix [dot(f[i].value,G[j])+f[i].divergence*H[j] for i in 1:3, j in 1:3]
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

        G = @SVector [j*T*g[i].value for i in 1:3]

        @SMatrix [dot(f[i].value, G[j]) for i in 1:3, j in 1:3]
end

kernel_in_bary(op::MWSingleLayer3D,
    test_local_space::RTRefSpace, trial_local_space::RTRefSpace,
    test_chart, trial_chart) = MWSL3DIntegrand(
        op, test_local_space, trial_local_space, test_chart, trial_chart)

kernel_in_bary(op::MWSingleLayer3D,
    test_local_space::RTRefSpace, trial_local_space::RTRefSpace,
    test_chart, trial_chart) = MWDL3DIntegrand(
        op, test_local_space, trial_local_space, test_chart, trial_chart)

function momintegrals!(op::MWSingleLayer3D,
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

    igd = MWSL3DIntegrand(test_triangular_element, trial_triangular_element,
        op, test_local_space, trial_local_space)
    G = sauterschwab_parameterized(igd, strat)
    for j ∈ 1:3, i ∈ 1:3
        out[i,j] += G[K[i],L[j]]
    end

    nothing
end






# function momintegrals!(op::MWDoubleLayer3D,
#     test_local_space::RTRefSpace, trial_local_space::RTRefSpace,
#     test_triangular_element, trial_triangular_element, out, strat::SauterSchwabStrategy)
#
#     I, J, K, L = SauterSchwabQuadrature.reorder(
#         test_triangular_element.vertices,
#         trial_triangular_element.vertices, strat)
#
#     test_triangular_element  = simplex(
#         test_triangular_element.vertices[I[1]],
#         test_triangular_element.vertices[I[2]],
#         test_triangular_element.vertices[I[3]])
#
#     trial_triangular_element = simplex(
#         trial_triangular_element.vertices[J[1]],
#         trial_triangular_element.vertices[J[2]],
#         trial_triangular_element.vertices[J[3]])
#
#     igd = MWDL3DIntegrand(
#         test_triangular_element, trial_triangular_element,
#         op, test_local_space, trial_local_space)
#     Q = sauterschwab_parameterized(igd, strat)
#     for j ∈ 1:3, i ∈ 1:3
#         out[i,j] += Q[K[i],L[j]]
#     end
# end
