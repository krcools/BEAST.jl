function pulled_back_integrand(op::HH3DSingleLayerFDBIO,
    test_local_space::LagrangeRefSpace{<:Any,0},
    trial_local_space::LagrangeRefSpace{<:Any,0},
    test_chart, trial_chart)

    (u,v) -> begin

        x = neighborhood(test_chart,u)
        y = neighborhood(trial_chart,v)

        f = test_local_space(x)
        g = trial_local_space(y)

        j = jacobian(x) * jacobian(y)

        α = alpha(op)
        γ = gamma(op)
        R = norm(cartesian(x)-cartesian(y))
        G = exp(-γ*R)/(4*π*R)

        αjG = α*G*j

        SA[f[1].value * αjG * g[1].value]
    end
end


function pulled_back_integrand(op::HH3DHyperSingularFDBIO,
    test_local_space::LagrangeRefSpace{<:Any,1},
    trial_local_space::LagrangeRefSpace{<:Any,1},
    test_chart, trial_chart)

    (u,v) -> begin

        x = neighborhood(test_chart,u)
        y = neighborhood(trial_chart,v)

        nx = normal(x)
        ny = normal(y)

        f = test_local_space(x)
        g = trial_local_space(y)

        j = jacobian(x) * jacobian(y)

        α = alpha(op)
        β = beta(op)
        γ = gamma(op)
        R = norm(cartesian(x)-cartesian(y))
        G = exp(-γ*R)/(4*π*R)

        αjG = ny*α*G*j
        βjG = β*G*j

        A = SA[ αjG*g[1].value, αjG*g[2].value, αjG*g[3].value]
        B = SA[ βjG*g[1].curl, βjG*g[2].curl, βjG*g[3].curl ] 

        SMatrix{3,3}((
            dot(nx*f[1].value, A[1]) + dot(f[1].curl, B[1]),
            dot(nx*f[2].value, A[1]) + dot(f[2].curl, B[1]),
            dot(nx*f[3].value, A[1]) + dot(f[3].curl, B[1]),
            dot(nx*f[1].value, A[2]) + dot(f[1].curl, B[2]),
            dot(nx*f[2].value, A[2]) + dot(f[2].curl, B[2]),
            dot(nx*f[3].value, A[2]) + dot(f[3].curl, B[2]),
            dot(nx*f[1].value, A[3]) + dot(f[1].curl, B[3]),
            dot(nx*f[2].value, A[3]) + dot(f[2].curl, B[3]),
            dot(nx*f[3].value, A[3]) + dot(f[3].curl, B[3])))
    end
end

function pulled_back_integrand(op::HH3DDoubleLayerFDBIO,
    test_local_space::LagrangeRefSpace{<:Any,0},
    trial_local_space::LagrangeRefSpace{<:Any,1},
    test_chart, trial_chart)

    (u,v) -> begin

        x = neighborhood(test_chart,u)
        y = neighborhood(trial_chart,v)

        ny = normal(y)
        f = test_local_space(x)
        g = trial_local_space(y)

        j = jacobian(x) * jacobian(y)

        α = alpha(op)
        γ = gamma(op)

        r = cartesian(x) - cartesian(y)
        R = norm(r)
        G = exp(-γ*R)/(4*π*R)
        inv_R = 1/R
        ∇G = -(γ + inv_R) * G * inv_R * r
        αnyj∇G = dot(ny,-α*∇G*j)

        SMatrix{1,3}(
            f[1].value * αnyj∇G * g[1].value,
            f[1].value * αnyj∇G * g[2].value,
            f[1].value * αnyj∇G * g[3].value)
    end
end

function pulled_back_integrand(op::HH3DDoubleLayerTransposedFDBIO,
    test_local_space::LagrangeRefSpace{<:Any,1},
    trial_local_space::LagrangeRefSpace{<:Any,0},
    test_chart, trial_chart)

    (u,v) -> begin

        x = neighborhood(test_chart,u)
        y = neighborhood(trial_chart,v)

        nx = normal(x)
        f = test_local_space(x)
        g = trial_local_space(y)

        j = jacobian(x) * jacobian(y)

        α = alpha(op)
        γ = gamma(op)

        r = cartesian(x) - cartesian(y)
        R = norm(r)
        G = exp(-γ*R)/(4*π*R)
        inv_R = 1/R
        ∇G = -(γ + inv_R) * G * inv_R * r
        αnxj∇G = dot(nx,α*∇G*j)

        SMatrix{3,1}(
            f[1].value * αnxj∇G * g[1].value,
            f[2].value * αnxj∇G * g[1].value,
            f[3].value * αnxj∇G * g[1].value)
    end
end

function momintegrals!(op::HH3DSingleLayerFDBIO,
    test_local_space::LagrangeRefSpace{<:Any,0}, trial_local_space::LagrangeRefSpace{<:Any,0},
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
    igd = pulled_back_integrand(op, test_local_space, trial_local_space,
        test_triangular_element, trial_triangular_element)
    G = SauterSchwabQuadrature.sauterschwab_parameterized(igd, strat)
    out[1,1] += G[1,1]

    nothing
end


function momintegrals!(op::HH3DHyperSingularFDBIO,
    test_local_space::LagrangeRefSpace{<:Any,1}, trial_local_space::LagrangeRefSpace{<:Any,1},
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

    test_sign = Combinatorics.levicivita(I)
    trial_sign = Combinatorics.levicivita(J)
    σ = test_sign * trial_sign

    igd = pulled_back_integrand(op, test_local_space, trial_local_space,
        test_triangular_element, trial_triangular_element)
    G = SauterSchwabQuadrature.sauterschwab_parameterized(igd, strat)
    for j ∈ 1:3, i ∈ 1:3
        out[i,j] += G[K[i],L[j]] * σ 
    end

    nothing
end

function momintegrals!(op::Helmholtz3DOp,
    test_local_space::LagrangeRefSpace{<:Any,0}, trial_local_space::LagrangeRefSpace{<:Any,1},
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

    trial_sign = Combinatorics.levicivita(J)

    σ = trial_sign

    igd = pulled_back_integrand(op, test_local_space, trial_local_space,
        test_triangular_element, trial_triangular_element)
    G = SauterSchwabQuadrature.sauterschwab_parameterized(igd, strat)

    for i ∈ 1:3
        out[1,i] += G[L[i]] * σ
    end

    nothing
end

function momintegrals!(op::Helmholtz3DOp,
    test_local_space::LagrangeRefSpace{<:Any,1}, trial_local_space::LagrangeRefSpace{<:Any,0},
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

        test_sign = Combinatorics.levicivita(I)
        σ = test_sign

    igd = pulled_back_integrand(op, test_local_space, trial_local_space,
        test_triangular_element, trial_triangular_element)
    G = SauterSchwabQuadrature.sauterschwab_parameterized(igd, strat)

    for i ∈ 1:3
        out[i,1] += G[K[i]] * σ
    end

    nothing
end
