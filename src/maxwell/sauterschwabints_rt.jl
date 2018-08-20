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

        jx = jacobian(x)
        jy = jacobian(y)
        j = jx*jy

        αjG = α*j*G
        βjG = β*j*G

        G1 = αjG*g[1][1]
        G2 = αjG*g[2][1]
        G3 = αjG*g[3][1]

        H1 = βjG*g[1][2]
        H2 = βjG*g[2][2]
        H3 = βjG*g[3][2]

        f1 = f[1][1]
        f2 = f[2][1]
        f3 = f[3][1]

        c1 = f[1][2]
        c2 = f[2][2]
        c3 = f[3][2]

        SMatrix{3,3}(
            dot(f1,G1) + c1*H1,
            dot(f2,G1) + c2*H1,
            dot(f3,G1) + c3*H1,
            dot(f1,G2) + c1*H2,
            dot(f2,G2) + c2*H2,
            dot(f3,G2) + c3*H2,
            dot(f1,G3) + c1*H3,
            dot(f2,G3) + c2*H3,
            dot(f3,G3) + c3*H3,)
end

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
    for j ∈ 1:3
        for i ∈ 1:3
            out[i,j] += G[K[i],L[j]]
        end
    end

    nothing
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

        GG = -(γ + 1/R) * G / R * r
        T = @SMatrix [
            0 -GG[3] GG[2]
            GG[3] 0 -GG[1]
            -GG[2] GG[1] 0 ]

        f = igd.test_local_space(x)
        g = igd.trial_local_space(y)

        jx = jacobian(x)
        jy = jacobian(y)
        j = jx*jy

        G1 = j*T*g[1][1]
        G2 = j*T*g[2][1]
        G3 = j*T*g[3][1]

        f1 = f[1][1]
        f2 = f[2][1]
        f3 = f[3][1]

        SMatrix{3,3}(
            dot(f1,G1),
            dot(f2,G1),
            dot(f3,G1),
            dot(f1,G2),
            dot(f2,G2),
            dot(f3,G2),
            dot(f1,G3),
            dot(f2,G3),
            dot(f3,G3),)
end


function momintegrals!(op::MWDoubleLayer3D,
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

    igd = MWDL3DIntegrand(test_triangular_element, trial_triangular_element,
        op, test_local_space, trial_local_space)
    Q = sauterschwab_parameterized(igd, strat)
    for j ∈ 1:3
        for i ∈ 1:3
            out[i,j] += Q[K[i],L[j]]
        end
    end
    nothing
end
