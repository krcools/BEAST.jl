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

function _integrands_gen(::Type{U}, ::Type{V}) where {U<:SVector{N}, V<:SVector{M}} where {M,N}
    ex = :(SMatrix{N,M}(()))
    for m in 1:M
        for n in 1:N
            # push!(ex.args[2].args, :(dot(a[$n], b[$m])))
            push!(ex.args[2].args, :(f(a[$n], b[$m])))
        end
    end
    return ex
end

@generated function _integrands(f, a::SVector{N}, b::SVector{M}) where {M,N}
    ex = _integrands_gen(a,b)
    # println(ex)
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
end
