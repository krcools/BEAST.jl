struct CurlSingleLayerDP3D{T,U} <: MaxwellOperator3D
    gamma::T
    alpha::U
end

function integrand(op::CurlSingleLayerDP3D, kernel_vals,
    test_vals, test_nbhd, trial_vals, trial_nbhd)

    α = op.alpha

    gx = test_vals.value
    fy = trial_vals.value

    ∇G = kernel_vals.gradgreen
    nx = normal(test_nbhd)

    return -α * dot(nx × gx, ∇G * fy)
end

function momintegrals!(op::CurlSingleLayerDP3D,
    test_local_space::RTRefSpace, trial_local_space::LagrangeRefSpace,
    test_triangular_element, trial_triangular_element, out, strat::SauterSchwabStrategy)

    I, J, K, L = SauterSchwabQuadrature.reorder(
        test_triangular_element.vertices,
        trial_triangular_element.vertices, strat)

    test_triangular_element  = simplex(test_triangular_element.vertices[I]...)
    trial_triangular_element = simplex(trial_triangular_element.vertices[J]...)

    function igd(u,v)
        α = op.alpha
        γ = op.gamma

        x = neighborhood(test_triangular_element,u)
        y = neighborhood(trial_triangular_element,v)

        nx = normal(x)

        r = cartesian(x) - cartesian(y)
        R = norm(r)
        G = exp(-γ*R)/(4π*R)
        GG = -(γ + 1/R) * G / R * r

        f = test_local_space(x)
        g = trial_local_space(y)

        jx = jacobian(x)
        jy = jacobian(y)
        j = jx*jy

        αjGG = α*j*GG

        R = @SVector[αjGG*g[i].value for i in 1:3]
        return @SMatrix[dot(f[i].value × nx, R[j]) for i in 1:3, j in 1:3]
    end

    Q = sauterschwab_parameterized(igd, strat)
    for j ∈ 1:3
        for i ∈ 1:3
            out[i,j] += Q[K[i],L[j]]
        end
    end

    nothing
end

function quadrule(op::BEAST.CurlSingleLayerDP3D,
    test_local_space::BEAST.RTRefSpace, trial_local_space::BEAST.RTRefSpace,
    test_index, test_chart, trial_index, trial_chart, quadrature_data)

    qrdf(op, test_local_space, trial_local_space, test_index, test_chart, trial_index, trial_chart, quadrature_data)
end
