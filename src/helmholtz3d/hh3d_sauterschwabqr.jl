function pulled_back_integrand(op::HH3DSingleLayerFDBIO,
    test_local_space::LagrangeRefSpace,
    trial_local_space::LagrangeRefSpace,
    test_chart, trial_chart)

    (u,v) -> begin

        x = neighborhood(test_chart,u)
        y = neighborhood(trial_chart,v)

        f = test_local_space(x)
        g = trial_local_space(y)

        j = jacobian(x) * jacobian(y)

        α = op.alpha
        γ = gamma(op)
        R = norm(cartesian(x)-cartesian(y))
        G = exp(-γ*R)/(4*π*R)

        αjG = α*G*j

        SMatrix{length(f),length(g)}((f[j].value * αjG * g[i].value for i in 1:length(g) for j in 1:length(f) )...)
    end
end


function pulled_back_integrand(op::HH3DHyperSingularFDBIO,
    test_local_space::LagrangeRefSpace,
    trial_local_space::LagrangeRefSpace,
    test_chart, trial_chart)

    (u,v) -> begin

        x = neighborhood(test_chart,u)
        y = neighborhood(trial_chart,v)

        nx = normal(x)
        ny = normal(y)

        f = test_local_space(x)
        g = trial_local_space(y)

        j = jacobian(x) * jacobian(y)

        α = op.alpha
        β = op.beta
        γ = gamma(op)
        R = norm(cartesian(x)-cartesian(y))
        G = exp(-γ*R)/(4*π*R)

        αjG = ny*α*G*j
        βjG = β*G*j

        A = SA[(αjG*g[i].value for i in 1:length(g))...]
        B = SA[(βjG*g[i].curl for i in 1:length(g))...]

        SMatrix{length(f),length(g)}((((dot(nx*f[j].value,A[i])+dot(f[j].curl,B[i])) for i in 1:length(g) for j in 1:length(f))...))
    end
end

function pulled_back_integrand(op::HH3DDoubleLayerFDBIO,
    test_local_space::LagrangeRefSpace,
    trial_local_space::LagrangeRefSpace,
    test_chart, trial_chart)

    (u,v) -> begin

        x = neighborhood(test_chart,u)
        y = neighborhood(trial_chart,v)

        ny = normal(y)
        f = test_local_space(x)
        g = trial_local_space(y)

        j = jacobian(x) * jacobian(y)

        α = op.alpha
        γ = gamma(op)

        r = cartesian(x) - cartesian(y)
        R = norm(r)
        G = exp(-γ*R)/(4*π*R)
        inv_R = 1/R
        ∇G = -(γ + inv_R) * G * inv_R * r
        αnyj∇G = dot(ny,-α*∇G*j)

        SMatrix{length(f),length(g)}((f[j].value * αnyj∇G * g[i].value  for i in 1:length(g) for j in 1:length(f))...)
    end
end

function pulled_back_integrand(op::HH3DDoubleLayerTransposedFDBIO,
    test_local_space::LagrangeRefSpace,
    trial_local_space::LagrangeRefSpace,
    test_chart, trial_chart)

    (u,v) -> begin

        x = neighborhood(test_chart,u)
        y = neighborhood(trial_chart,v)

        nx = normal(x)
        f = test_local_space(x)
        g = trial_local_space(y)

        j = jacobian(x) * jacobian(y)

        α = op.alpha
        γ = gamma(op)

        r = cartesian(x) - cartesian(y)
        R = norm(r)
        G = exp(-γ*R)/(4*π*R)
        inv_R = 1/R
        ∇G = -(γ + inv_R) * G * inv_R * r
        αnxj∇G = dot(nx,α*∇G*j)

        SMatrix{length(f),length(g)}((f[j].value * αnxj∇G * g[i].value for i in 1:length(g) for j in 1:length(f))...)
    end
end

function momintegrals!(op::Helmholtz3DOp,
    test_local_space::LagrangeRefSpace{<:Any,0}, trial_local_space::LagrangeRefSpace{<:Any,0},
    test_triangular_element, trial_triangular_element, out, strat::SauterSchwabStrategy)

    out_old = copy(out)
    momintegrals!_old(op,test_local_space,trial_local_space,test_triangular_element,trial_triangular_element,out_old,strat)
    I, J, K, L = SauterSchwabQuadrature.reorder(
        test_triangular_element.vertices,
        trial_triangular_element.vertices, strat)

    # test_triangular_elementex  = simplex(
    #     test_triangular_element.vertices[I[1]],
    #     test_triangular_element.vertices[I[2]],
    #     test_triangular_element.vertices[I[3]])

    # trial_triangular_elementex = simplex(
    #     trial_triangular_element.vertices[J[1]],
    #     trial_triangular_element.vertices[J[2]],
    #     trial_triangular_element.vertices[J[3]])

    #     test_sign = Combinatorics.levicivita(I)
    #     trial_sign = Combinatorics.levicivita(J)
    #     σ = momintegrals_sign(op, test_sign, trial_sign)
    #     @assert sign(dot(normal(test_triangular_element),normal(test_triangular_elementex))) == test_sign
    # test_triangular_element = CompScienceMeshes.Simplex(test_triangular_elementex.vertices,
    #                                                     test_triangular_elementex.tangents,
    #                                                     test_triangular_element.normals,
    #                                                     test_triangular_elementex.volume)
    # trial_triangular_element = CompScienceMeshes.Simplex(trial_triangular_elementex.vertices,
    #                                                 trial_triangular_elementex.tangents,
    #                                                 trial_triangular_element.normals,
    #                                                 trial_triangular_elementex.volume)


    igd = pulled_back_integrand(op, test_local_space, trial_local_space,
        test_triangular_element, trial_triangular_element)

    f = (u,v) -> igd(permute_barycentric(K,u),permute_barycentric(L,v))
    G = SauterSchwabQuadrature.sauterschwab_parameterized(f, strat)
    out[1,1] += G[1,1]# * σ

    @assert out_old ≈ out rtol = eps()
    nothing
end


function momintegrals!(op::Helmholtz3DOp,
    test_local_space::LagrangeRefSpace, trial_local_space::LagrangeRefSpace,
    test_triangular_element, trial_triangular_element, out, strat::SauterSchwabStrategy)
    out_old = copy(out)
    momintegrals!_old(op,test_local_space,trial_local_space,test_triangular_element,trial_triangular_element,out_old,strat)
    I, J, K, L = SauterSchwabQuadrature.reorder(
        test_triangular_element.vertices,
        trial_triangular_element.vertices, strat)

    # test_triangular_elementex  = simplex(
    #     test_triangular_element.vertices[I[1]],
    #     test_triangular_element.vertices[I[2]],
    #     test_triangular_element.vertices[I[3]])

    # trial_triangular_elementex = simplex(
    #     trial_triangular_element.vertices[J[1]],
    #     trial_triangular_element.vertices[J[2]],
    #     trial_triangular_element.vertices[J[3]])

    #     test_sign = Combinatorics.levicivita(I)
    #     trial_sign = Combinatorics.levicivita(J)
    #     σ = momintegrals_sign(op, test_sign, trial_sign)
    #     @assert sign(dot(normal(test_triangular_element),normal(test_triangular_elementex))) == test_sign
    # test_triangular_element = CompScienceMeshes.Simplex(test_triangular_elementex.vertices,
    #                                                     test_triangular_elementex.tangents,
    #                                                     test_triangular_element.normals,
    #                                                     test_triangular_elementex.volume)
    # trial_triangular_element = CompScienceMeshes.Simplex(trial_triangular_elementex.vertices,
    #                                                 trial_triangular_elementex.tangents,
    #                                                 trial_triangular_element.normals,
    #                                                 trial_triangular_elementex.volume)

    igd = pulled_back_integrand(op, test_local_space, trial_local_space,
        test_triangular_element, trial_triangular_element)
    f = (u,v) -> igd(permute_barycentric(K,u),permute_barycentric(L,v))
    G = SauterSchwabQuadrature.sauterschwab_parameterized(f, strat)
    for j ∈ 1:3, i ∈ 1:3
       # out[i,j] += G[K[i],L[j]]# * σ
        out[i,j] += G[i,j]# * σ
    end
    @assert out_old ≈ out rtol = eps()
    nothing
end

function momintegrals!(op::Helmholtz3DOp,
    test_local_space::LagrangeRefSpace{<:Any,0}, trial_local_space::LagrangeRefSpace{<:Any,1},
    test_triangular_element, trial_triangular_element, out, strat::SauterSchwabStrategy)
    out_old = copy(out)
    momintegrals!_old(op,test_local_space,trial_local_space,test_triangular_element,trial_triangular_element,out_old,strat)
    I, J, K, L = SauterSchwabQuadrature.reorder(
        test_triangular_element.vertices,
        trial_triangular_element.vertices, strat)
    
    # test_triangular_elementex  = simplex(
    #     test_triangular_element.vertices[I[1]],
    #     test_triangular_element.vertices[I[2]],
    #     test_triangular_element.vertices[I[3]])

    # trial_triangular_elementex = simplex(
    #     trial_triangular_element.vertices[J[1]],
    #     trial_triangular_element.vertices[J[2]],
    #     trial_triangular_element.vertices[J[3]])

    #     test_sign = Combinatorics.levicivita(I)
    #     trial_sign = Combinatorics.levicivita(J)
    #     σ = momintegrals_sign(op, test_sign, trial_sign)
    #     @assert sign(dot(normal(test_triangular_element),normal(test_triangular_elementex))) == test_sign
    # test_triangular_element = CompScienceMeshes.Simplex(test_triangular_elementex.vertices,
    #                                                     test_triangular_elementex.tangents,
    #                                                     test_triangular_element.normals,
    #                                                     test_triangular_elementex.volume)
    # trial_triangular_element = CompScienceMeshes.Simplex(trial_triangular_elementex.vertices,
    #                                                 trial_triangular_elementex.tangents,
    #                                                 trial_triangular_element.normals,
    #                                                 trial_triangular_elementex.volume)

    igd = pulled_back_integrand(op, test_local_space, trial_local_space,
        test_triangular_element, trial_triangular_element)
    f = (u,v) -> igd(permute_barycentric(K,u),permute_barycentric(L,v))    
    G = SauterSchwabQuadrature.sauterschwab_parameterized(f, strat)

    for i ∈ 1:3
        #out[1,i] += G[L[i]]# * σ
        out[1,i] += G[i]# * σ
    end
    @assert out_old ≈ out rtol = eps()
    nothing
end

function momintegrals!(op::Helmholtz3DOp,
    test_local_space::LagrangeRefSpace{<:Any,1}, trial_local_space::LagrangeRefSpace{<:Any,0},
    test_triangular_element, trial_triangular_element, out, strat::SauterSchwabStrategy)
    out_old = copy(out)
    momintegrals!_old(op,test_local_space,trial_local_space,test_triangular_element,trial_triangular_element,out_old,strat)
    I, J, K, L = SauterSchwabQuadrature.reorder(
        test_triangular_element.vertices,
        trial_triangular_element.vertices, strat)

    # test_triangular_element  = simplex(
    #     test_triangular_element.vertices[I[1]],
    #     test_triangular_element.vertices[I[2]],
    #     test_triangular_element.vertices[I[3]])

    # trial_triangular_element = simplex(
    #     trial_triangular_element.vertices[J[1]],
    #     trial_triangular_element.vertices[J[2]],
    #     trial_triangular_element.vertices[J[3]])
    # # test_triangular_element = CompScienceMeshes.Simplex(test_triangular_element1.vertices,
    # #                                                     test_triangular_element1.tangents,
    # #                                                     test_triangular_element.normals,
    # #                                                     test_triangular_element1.volume)
    # # trial_triangular_element = CompScienceMeshes.Simplex(trial_triangular_element1.vertices,
    # #                                                     trial_triangular_element1.tangents,
    # #                                                     trial_triangular_element.normals,
    # #                                                     trial_triangular_element1.volume)
    #     test_sign = Combinatorics.levicivita(I)
    #     trial_sign = Combinatorics.levicivita(J)

    #     σ = momintegrals_sign(op, test_sign, trial_sign)

    igd = pulled_back_integrand(op, test_local_space, trial_local_space,
        test_triangular_element, trial_triangular_element)
    f = (u,v) -> igd(permute_barycentric(K,u),permute_barycentric(L,v))   
    G = SauterSchwabQuadrature.sauterschwab_parameterized(f, strat)

    for i ∈ 1:3
        #out[i,1] += G[K[i]] * σ
        out[i,1] += G[i] 
    end
    @assert out_old ≈ out rtol = eps()
    nothing
end

function momintegrals_sign(op::HH3DSingleLayerFDBIO, test_sign, trial_sign)
    return 1
end
function momintegrals_sign(op::HH3DDoubleLayerFDBIO, test_sign, trial_sign)
    return trial_sign
end
function momintegrals_sign(op::HH3DDoubleLayerTransposedFDBIO, test_sign, trial_sign)
    return test_sign
end
function momintegrals_sign(op::HH3DHyperSingularFDBIO, test_sign, trial_sign)
    return test_sign * trial_sign
end



function momintegrals!_old(op::Helmholtz3DOp,
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

    test_sign = Combinatorics.levicivita(I)
    trial_sign = Combinatorics.levicivita(J)
    σ = momintegrals_sign(op, test_sign, trial_sign)

    igd = pulled_back_integrand(op, test_local_space, trial_local_space,
        test_triangular_element, trial_triangular_element)
    G = SauterSchwabQuadrature.sauterschwab_parameterized(igd, strat)
    out[1,1] += G[1,1] * σ

    nothing
end


function momintegrals!_old(op::Helmholtz3DOp,
    test_local_space::LagrangeRefSpace, trial_local_space::LagrangeRefSpace,
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

    σ = momintegrals_sign(op, test_sign, trial_sign)

    igd = pulled_back_integrand(op, test_local_space, trial_local_space,
        test_triangular_element, trial_triangular_element)
    G = SauterSchwabQuadrature.sauterschwab_parameterized(igd, strat)
    for j ∈ 1:3, i ∈ 1:3
        out[i,j] += G[K[i],L[j]] * σ
    end

    nothing
end

function momintegrals!_old(op::Helmholtz3DOp,
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

    test_sign = Combinatorics.levicivita(I)
    trial_sign = Combinatorics.levicivita(J)

    σ = momintegrals_sign(op, test_sign, trial_sign)

    igd = pulled_back_integrand(op, test_local_space, trial_local_space,
        test_triangular_element, trial_triangular_element)
    G = SauterSchwabQuadrature.sauterschwab_parameterized(igd, strat)

    for i ∈ 1:3
        out[1,i] += G[L[i]] * σ
    end

    nothing
end

function momintegrals!_old(op::Helmholtz3DOp,
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
        trial_sign = Combinatorics.levicivita(J)

        σ = momintegrals_sign(op, test_sign, trial_sign)

    igd = pulled_back_integrand(op, test_local_space, trial_local_space,
        test_triangular_element, trial_triangular_element)
    G = SauterSchwabQuadrature.sauterschwab_parameterized(igd, strat)

    for i ∈ 1:3
        out[i,1] += G[K[i]] * σ
    end

    nothing
end