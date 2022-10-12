


const i4pi = 1 / (4pi)
function (igd::Integrand)(u,v)
    
    x = neighborhood(igd.test_chart,u)
    y = neighborhood(igd.trial_chart,v)
    
    f = igd.local_test_space(x)
    g = igd.local_trial_space(y)

    return jacobian(x) * jacobian(y) * igd(x,y,f,g)
end


function _integrands_leg_gen(f::Type{U}, g::Type{V}) where {U<:SVector{N}, V<:SVector{M}} where {M,N}
    ex = :(SMatrix{N,M}(()))
    for m in 1:M
        for n in 1:N
            push!(ex.args[2].args, :(integrand(op, kervals, f[$n], x, g[$m], y)))
        end
    end
    return ex
end

@generated _integrands_leg(op, kervals, f::SVector{N}, x, g::SVector{M}, y) where {M,N} = _integrands_leg_gen(f, g)

# Support for legacy kernels
function (igd::Integrand)(x,y,f,g)

    op = igd.operator
    kervals = kernelvals(op, x, y)
    _integrands_leg(op, kervals, f, x, g, y)

end


function (igd::Integrand{<:MWSingleLayer3D})(x,y,f,g)
    α = igd.operator.α
    β = igd.operator.β
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1 / R
    green = exp(-γ*R)*(i4pi*iR)

    αG = α * green
    βG = β * green

    _integrands(f,g) do fi,gj
        αG * dot(fi.value, gj.value) + βG * dot(fi.divergence, gj.divergence)
    end
end

function (igd::Integrand{<:MWSingleLayer3DReg})(x,y,f,g)
    α = igd.operator.α
    β = igd.operator.β
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    γR = γ*R
    # iR = 1 / R
    green = (expm1(-γR) + γR - 0.5*γR^2) / (4pi*R)

    αG = α * green
    βG = β * green

    _integrands(f,g) do fi,gj
        αG * dot(fi.value, gj.value) + βG * dot(fi.divergence, gj.divergence)
    end
end


function (igd::Integrand{<:MWDoubleLayer3D})(x,y,f,g)
    
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)
    gradgreen = -(γ + iR) * green * (iR * r)

    fvalue = getvalue(f)
    gvalue = getvalue(g)
    G = cross.(Ref(gradgreen), gvalue)
    return _krondot(fvalue, G)
end


function (igd::Integrand{<:MWDoubleLayer3DReg})(x,y,f,g)
    
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    γR = γ*R
    iR = 1/R
    expo = exp(-γR)
    green = (expo - 1 + γR - 0.5*γR^2) * (i4pi*iR)
    gradgreen = ( -(γR + 1)*expo + (1 - 0.5*γR^2) ) * (i4pi*iR^3) * r

    fvalue = getvalue(f)
    gvalue = getvalue(g)
    G = cross.(Ref(gradgreen), gvalue)
    return _krondot(fvalue, G)
end

function momintegrals_test_refines_trial!(out, op,
    test_functions, test_cell, test_chart,
    trial_functions, trial_cell, trial_chart,
    quadrule, quadstrat)

    test_local_space = refspace(test_functions)
    trial_local_space = refspace(trial_functions)

    momintegrals!(op, test_local_space, trial_local_space,
        test_chart, trial_chart, out, quadrule)
end

# const MWOperator3D = Union{MWSingleLayer3D, MWDoubleLayer3D}
function momintegrals_test_refines_trial!(out, op,
    test_functions, test_cell, test_chart,
    trial_functions, trial_cell, trial_chart,
    qr::SauterSchwabStrategy, quadstrat)

    test_local_space = refspace(test_functions)
    trial_local_space = refspace(trial_functions)

    test_mesh = geometry(test_functions)
    trial_mesh = geometry(trial_functions)

    parent_mesh = CompScienceMeshes.parent(test_mesh)
    trial_charts = [chart(test_mesh, p) for p in CompScienceMeshes.children(parent_mesh, trial_cell)]

    qd = quaddata(op, test_local_space, trial_local_space,
        [test_chart], trial_charts, quadstrat)

    for (q,chart) in enumerate(trial_charts)
        qr = quadrule(op, test_local_space, trial_local_space,
            1, test_chart, q ,chart, qd, quadstrat)

        Q = restrict(trial_local_space, trial_chart, chart)
        zlocal = zero(out)
        momintegrals!(op, test_local_space, trial_local_space,
            test_chart, chart, zlocal, qr)

        for j in 1:3
            for i in 1:3
                for k in 1:3
                    out[i,j] += zlocal[i,k] * Q[j,k]
end end end end end



function momintegrals_trial_refines_test!(out, op,
    test_functions, test_cell, test_chart,
    trial_functions, trial_cell, trial_chart,
    quadrule, quadstrat)

    test_local_space = refspace(test_functions)
    trial_local_space = refspace(trial_functions)

    momintegrals!(op, test_local_space, trial_local_space,
        test_chart, trial_chart, out, quadrule)
end


function momintegrals_trial_refines_test!(out, op,
    test_functions, test_cell, test_chart,
    trial_functions, trial_cell, trial_chart,
    qr::SauterSchwabStrategy, quadstrat)

    test_local_space = refspace(test_functions)
    trial_local_space = refspace(trial_functions)

    test_mesh = geometry(test_functions)
    trial_mesh = geometry(trial_functions)

    parent_mesh = CompScienceMeshes.parent(trial_mesh)
    test_charts = [chart(trial_mesh, p) for p in CompScienceMeshes.children(parent_mesh, test_cell)] 

    qd = quaddata(op, test_local_space, trial_local_space,
        test_charts, [trial_chart], quadstrat)

    for (p,chart) in enumerate(test_charts)
        qr = quadrule(op, test_local_space, trial_local_space,
            p, chart, 1, trial_chart, qd, quadstrat)

        Q = restrict(test_local_space, test_chart, chart)
        zlocal = zero(out)
        momintegrals!(op, test_local_space, trial_local_space,
            chart, trial_chart, zlocal, qr)

        for j in 1:3
            for i in 1:3
                for k in 1:3
                    out[i,j] += Q[i,k] * zlocal[k,j]
        end end end
    end
end