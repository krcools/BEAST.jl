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

        j = jacobian(x) * jacobian(y)

        αjG = α*j*G
        βjG = β*j*G

        G = @SVector[αjG*g[i].value for i in 1:6]
        H = @SVector[βjG*g[i].divergence for i in 1:6]

        @SMatrix[dot(f[i].value,G[j])+f[i].divergence*H[j] for i in 1:6, j in 1:6]
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
