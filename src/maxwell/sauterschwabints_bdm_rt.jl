

struct MWDL3DIntegrand3{C,O,L,M}
        test_triangular_element::C
        trial_triangular_element::C
        op::O
        test_local_space::L
        trial_local_space::M
end

function (igd::MWDL3DIntegrand3)(u,v)

        γ = igd.op.gamma

        x = neighborhood(igd.test_triangular_element,u)
        y = neighborhood(igd.trial_triangular_element,v)

        r = cartesian(x) - cartesian(y)
        R = norm(r)
        G = exp(-γ*R)/(4π*R)

        GG = -(γ + 1/R) * G / R * r
  
        f = igd.test_local_space(x)  #BDM
        g = igd.trial_local_space(y) #RT

        j = jacobian(x) * jacobian(y)

        G = @SVector [cross(GG,(j*g[1].value)), cross(GG,(j*g[2].value)), cross(GG,(j*g[3].value))]
       

        f1 = f[1].value
        f2 = f[2].value
        f3 = f[3].value
        f4 = f[4].value
        f5 = f[5].value
        f6 = f[6].value

        SMatrix{6,3}(
            dot(f1,G[1]),
            dot(f2,G[1]),
            dot(f3,G[1]),
            dot(f4,G[1]),
            dot(f5,G[1]),
            dot(f6,G[1]),
            dot(f1,G[2]),
            dot(f2,G[2]),
            dot(f3,G[2]),
            dot(f4,G[2]),
            dot(f5,G[2]),
            dot(f6,G[2]),
            dot(f1,G[3]),
            dot(f2,G[3]),
            dot(f3,G[3]),
            dot(f4,G[3]),
            dot(f5,G[3]),
            dot(f6,G[3]))
end


function momintegrals!(op::MWDoubleLayer3D,
    test_local_space::BDMRefSpace, trial_local_space::RTRefSpace,
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

    igd = MWDL3DIntegrand3(test_triangular_element, trial_triangular_element,
        op, test_local_space, trial_local_space)
    Q = SauterSchwabQuadrature.sauterschwab_parameterized(igd, strat)

    A = levicivita(K) == 1 ? @SVector[1,2] : @SVector[2,1]

    for j ∈ 1:3
        for i ∈ 1:3
            p = K[i]
            for m in 1:2
                a = A[m]
                out[2*(i-1)+m,j] += Q[2*(p-1)+a,L[j]]
            end
        end
    end
    nothing
end
