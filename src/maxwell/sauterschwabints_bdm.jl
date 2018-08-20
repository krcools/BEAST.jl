struct MWSL3DIntegrand2{C,O,L}
        test_triangular_element::C
        trial_triangular_element::C
        op::O
        test_local_space::L
        trial_local_space::L
end

function (igd::MWSL3DIntegrand2)(u,v)
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

        # G1 = αjG*g[1].value
        # G2 = αjG*g[2].value
        # G3 = αjG*g[3].value
        # G4 = αjG*g[4].value
        # G5 = αjG*g[5].value
        # G6 = αjG*g[6].value

        G = @SVector[αjG*g[i].value for i in 1:6]

        # H1 = βjG*g[1].divergence
        # H2 = βjG*g[2].divergence
        # H3 = βjG*g[3].divergence
        # H4 = βjG*g[4].divergence
        # H5 = βjG*g[5].divergence
        # H6 = βjG*g[6].divergence

        H = @SVector[βjG*g[i].divergence for i in 1:6]

        # f1 = f[1].value
        # f2 = f[2].value
        # f3 = f[3].value
        # f4 = f[4].value
        # f5 = f[5].value
        # f6 = f[6].value

        # c1 = f[1].divergence
        # c2 = f[2].divergence
        # c3 = f[3].divergence
        # c4 = f[4].divergence
        # c5 = f[5].divergence
        # c6 = f[6].divergence

        @SMatrix[dot(f[i].value,G[j])+f[i].divergence*H[j] for i in 1:6, j in 1:6]

        # SMatrix{6,6}(
        #     dot(f1,G1) + c1*H1,
        #     dot(f2,G1) + c2*H1,
        #     dot(f3,G1) + c3*H1,
        #     dot(f1,G2) + c1*H2,
        #     dot(f2,G2) + c2*H2,
        #     dot(f3,G2) + c3*H2,
        #     dot(f1,G3) + c1*H3,
        #     dot(f2,G3) + c2*H3,
        #     dot(f3,G3) + c3*H3,)
end

function momintegrals!(op::MWSingleLayer3D,
    test_local_space::BDMRefSpace, trial_local_space::BDMRefSpace,
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

    igd = MWSL3DIntegrand2(test_triangular_element, trial_triangular_element,
        op, test_local_space, trial_local_space)
    G = sauterschwab_parameterized(igd, strat)

    A = levicivita(K) == 1 ? @SVector[1,2] : @SVector[2,1]
    B = levicivita(L) == 1 ? @SVector[1,2] : @SVector[2,1]
    for j ∈ 1:3
        q = L[j]
        for i ∈ 1:3
            p = K[i]
            for n in 1:2
                b = B[n]
                for m in 1:2
                    a = A[m]
                    out[2*(i-1)+m,2*(j-1)+n] += G[2*(p-1)+a,2*(q-1)+b]
                end
            end
        end
    end

    nothing
end


struct MWDL3DIntegrand2{C,O,L}
        test_triangular_element::C
        trial_triangular_element::C
        op::O
        test_local_space::L
        trial_local_space::L
end

function (igd::MWDL3DIntegrand2)(u,v)

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
    test_local_space::BDMRefSpace, trial_local_space::BDMRefSpace,
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

    igd = MWDL3DIntegrand2(test_triangular_element, trial_triangular_element,
        op, test_local_space, trial_local_space)
    Q = sauterschwab_parameterized(igd, strat)
    for j ∈ 1:3
        for i ∈ 1:3
            out[i,j] += Q[K[i],L[j]]
        end
    end
    nothing
end
