struct Integrand{Op,LSt,LSb,Elt,Elb}
    operator::Op
    local_test_space::LSt
    local_trial_space::LSb
    test_chart::Elt
    trial_chart::Elb
end

getvalue(a::SVector{N}) where {N} = SVector{N}(getvalue(a.data))
getvalue(a::NTuple{1}) = (a[1].value,)
getvalue(a::NTuple{N}) where {N} = tuple(a[1].value, getvalue(Base.tail(a))...)

getdivergence(a::SVector{N}) where {N} = SVector{N}(getdivergence(a.data))
getdivergence(a::NTuple{1}) = (a[1].divergence,)
getdivergence(a::NTuple{N}) where {N} = tuple(a[1].divergence, getdivergence(Base.tail(a))...)

function _krondot_gen(a::Type{U}, b::Type{V}) where {U<:SVector{N}, V<:SVector{M}} where {M,N}
    ex = :(SMatrix{N,M}(()))
    for m in 1:M
        for n in 1:N
            push!(ex.args[2].args, :(dot(a[$n], b[$m])))
        end
    end
    return ex
end

@generated function _krondot(a::SVector{N}, b::SVector{M}) where {M,N}
    ex = _krondot_gen(a,b)
    return ex
end

function momintegrals!(op::Operator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_chart, trial_chart, out, strat::SauterSchwabStrategy)

    I, J, _, _ = SauterSchwabQuadrature.reorder(
        test_chart.vertices,
        trial_chart.vertices, strat)

    test_chart  = simplex(
        test_chart.vertices[I[1]],
        test_chart.vertices[I[2]],
        test_chart.vertices[I[3]])

    trial_chart = simplex(
        trial_chart.vertices[J[1]],
        trial_chart.vertices[J[2]],
        trial_chart.vertices[J[3]])

    igd = Integrand(op, test_local_space, trial_local_space, test_chart, trial_chart)
    G = SauterSchwabQuadrature.sauterschwab_parameterized(igd, strat)

    K = dof_permutation(test_local_space, I)
    L = dof_permutation(trial_local_space, J)
    for i in 1:numfunctions(test_local_space)
        for j in 1:numfunctions(trial_local_space)
            out[i,j] = G[K[i],L[j]]
    end end

    nothing
end


function (igd::Integrand{<:DoubleLayerRotatedMW3D})(u,v)

    x = neighborhood(igd.test_chart,u)
    y = neighborhood(igd.trial_chart,v)
    j = jacobian(x) * jacobian(y)
    nx = normal(x)

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    γ = igd.operator.gamma
    G = exp(-γ*R)/(4π*R)
    K = -(γ + iR) * G * (iR * r)

    f = igd.local_test_space(x)
    g = igd.local_trial_space(y)

    fvalue = getvalue(f)
    gvalue = getvalue(g)

    jKg = cross.(Ref(K), j*gvalue)
    jnxKg = cross.(Ref(nx), jKg)
    return _krondot(fvalue, jnxKg)

    # Kg = [cross(K, j*gi.value) for gi in g]
    # return [dot(fj.value, cross(nx, Kgi)) for fj in f, Kgi in Kg]
end


# struct MWSL3DIntegrand2{C,O,L}
#         test_triangular_element::C
#         trial_triangular_element::C
#         op::O
#         test_local_space::L
#         trial_local_space::L
# end

# function (igd::MWSL3DIntegrand2)(u,v)
#         α = igd.op.α
#         β = igd.op.β
#         γ = igd.op.gamma

#         x = neighborhood(igd.test_triangular_element,u)
#         y = neighborhood(igd.trial_triangular_element,v)

#         r = cartesian(x) - cartesian(y)
#         R = norm(r)
#         G = exp(-γ*R)/(4π*R)

#         f = igd.test_local_space(x)
#         g = igd.trial_local_space(y)

#         j = jacobian(x) * jacobian(y)

#         αjG = α*j*G
#         βjG = β*j*G

#         G = @SVector[αjG*g[i].value for i in 1:6]
#         H = @SVector[βjG*g[i].divergence for i in 1:6]

#         SMatrix{6,6}(tuple(
#             dot(f[1].value,G[1])+f[1].divergence*H[1],
#             dot(f[2].value,G[1])+f[2].divergence*H[1],
#             dot(f[3].value,G[1])+f[3].divergence*H[1],
#             dot(f[4].value,G[1])+f[4].divergence*H[1],
#             dot(f[5].value,G[1])+f[5].divergence*H[1],
#             dot(f[6].value,G[1])+f[6].divergence*H[1],
#             dot(f[1].value,G[2])+f[1].divergence*H[2],
#             dot(f[2].value,G[2])+f[2].divergence*H[2],
#             dot(f[3].value,G[2])+f[3].divergence*H[2],
#             dot(f[4].value,G[2])+f[4].divergence*H[2],
#             dot(f[5].value,G[2])+f[5].divergence*H[2],
#             dot(f[6].value,G[2])+f[6].divergence*H[2],
#             dot(f[1].value,G[3])+f[1].divergence*H[3],
#             dot(f[2].value,G[3])+f[2].divergence*H[3],
#             dot(f[3].value,G[3])+f[3].divergence*H[3],
#             dot(f[4].value,G[3])+f[4].divergence*H[3],
#             dot(f[5].value,G[3])+f[5].divergence*H[3],
#             dot(f[6].value,G[3])+f[6].divergence*H[3],
#             dot(f[1].value,G[4])+f[1].divergence*H[4],
#             dot(f[2].value,G[4])+f[2].divergence*H[4],
#             dot(f[3].value,G[4])+f[3].divergence*H[4],
#             dot(f[4].value,G[4])+f[4].divergence*H[4],
#             dot(f[5].value,G[4])+f[5].divergence*H[4],
#             dot(f[6].value,G[4])+f[6].divergence*H[4],
#             dot(f[1].value,G[5])+f[1].divergence*H[5],
#             dot(f[2].value,G[5])+f[2].divergence*H[5],
#             dot(f[3].value,G[5])+f[3].divergence*H[5],
#             dot(f[4].value,G[5])+f[4].divergence*H[5],
#             dot(f[5].value,G[5])+f[5].divergence*H[5],
#             dot(f[6].value,G[5])+f[6].divergence*H[5],
#             dot(f[1].value,G[6])+f[1].divergence*H[6],
#             dot(f[2].value,G[6])+f[2].divergence*H[6],
#             dot(f[3].value,G[6])+f[3].divergence*H[6],
#             dot(f[4].value,G[6])+f[4].divergence*H[6],
#             dot(f[5].value,G[6])+f[5].divergence*H[6],
#             dot(f[6].value,G[6])+f[6].divergence*H[6],
#         ))
# end

# function momintegrals!(op::MWSingleLayer3D,
#     test_local_space::BDMRefSpace, trial_local_space::BDMRefSpace,
#     test_triangular_element, trial_triangular_element, out, strat::SauterSchwabStrategy)

#     I, J, _, _ = SauterSchwabQuadrature.reorder(
#         test_triangular_element.vertices,
#         trial_triangular_element.vertices, strat)

#     test_triangular_element  = simplex(
#         test_triangular_element.vertices[I[1]],
#         test_triangular_element.vertices[I[2]],
#         test_triangular_element.vertices[I[3]])

#     trial_triangular_element = simplex(
#         trial_triangular_element.vertices[J[1]],
#         trial_triangular_element.vertices[J[2]],
#         trial_triangular_element.vertices[J[3]])

#     igd = MWSL3DIntegrand2(test_triangular_element, trial_triangular_element,
#         op, test_local_space, trial_local_space)
#     G = SauterSchwabQuadrature.sauterschwab_parameterized(igd, strat)

#     K = dof_permutation(test_local_space, I)
#     L = dof_permutation(trial_local_space, J)
#     for i in 1:numfunctions(test_local_space)
#         for j in 1:numfunctions(trial_local_space)
#             out[i,j] = G[K[i],L[j]]
#     end end

#     nothing
# end


# struct MWDL3DIntegrand2{C,O,L}
#         test_triangular_element::C
#         trial_triangular_element::C
#         op::O
#         test_local_space::L
#         trial_local_space::L
# end

# function (igd::MWDL3DIntegrand2)(u,v)

#         γ = igd.op.gamma

#         x = neighborhood(igd.test_triangular_element,u)
#         y = neighborhood(igd.trial_triangular_element,v)

#         r = cartesian(x) - cartesian(y)
#         R = norm(r)
#         G = exp(-γ*R)/(4π*R)

#         GG = -(γ + 1/R) * G / R * r
#         T = @SMatrix [
#             0 -GG[3] GG[2]
#             GG[3] 0 -GG[1]
#             -GG[2] GG[1] 0 ]

#         f = igd.test_local_space(x)
#         g = igd.trial_local_space(y)

#         jx = jacobian(x)
#         jy = jacobian(y)
#         j = jx*jy

#         G1 = j*T*g[1][1]
#         G2 = j*T*g[2][1]
#         G3 = j*T*g[3][1]

#         f1 = f[1][1]
#         f2 = f[2][1]
#         f3 = f[3][1]

#         SMatrix{3,3}(
#             dot(f1,G1),
#             dot(f2,G1),
#             dot(f3,G1),
#             dot(f1,G2),
#             dot(f2,G2),
#             dot(f3,G2),
#             dot(f1,G3),
#             dot(f2,G3),
#             dot(f3,G3),)
# end


# function momintegrals!(op::MWDoubleLayer3D,
#     test_local_space::BDMRefSpace, trial_local_space::BDMRefSpace,
#     test_triangular_element, trial_triangular_element, out, strat::SauterSchwabStrategy)

#     I, J, K, L = SauterSchwabQuadrature.reorder(
#         test_triangular_element.vertices,
#         trial_triangular_element.vertices, strat)

#     test_triangular_element  = simplex(
#         test_triangular_element.vertices[I[1]],
#         test_triangular_element.vertices[I[2]],
#         test_triangular_element.vertices[I[3]])

#     trial_triangular_element = simplex(
#         trial_triangular_element.vertices[J[1]],
#         trial_triangular_element.vertices[J[2]],
#         trial_triangular_element.vertices[J[3]])

#     igd = MWDL3DIntegrand2(test_triangular_element, trial_triangular_element,
#         op, test_local_space, trial_local_space)
#     Q = sauterschwab_parameterized(igd, strat)
#     for j ∈ 1:3
#         for i ∈ 1:3
#             out[i,j] += Q[K[i],L[j]]
#         end
#     end
#     nothing
# end