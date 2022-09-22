


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

const MWOperator3D = Union{MWSingleLayer3D, MWDoubleLayer3D}
function momintegrals_nested!(op::MWOperator3D,
    test_local_space::RTRefSpace, trial_local_space::RTRefSpace,
    test_chart, trial_chart, out, strat::SauterSchwabStrategy, quadstrat)

    # 1. Refine the trial_chart
    p1, p2, p3 = trial_chart.vertices

    # TODO: generalise this to include more general refinements
    e1 = cartesian(neighborhood(trial_chart, (0,1/2)))
    e2 = cartesian(neighborhood(trial_chart, (1/2,0)))
    e3 = cartesian(neighborhood(trial_chart, (1/2,1/2)))

    ct = cartesian(center(trial_chart))

    refined_trial_chart = SA[
        simplex(ct, p1, e3),
        simplex(ct, e3, p2),
        simplex(ct, p2, e1),
        simplex(ct, e1, p3),
        simplex(ct, p3, e2),
        simplex(ct, e2, p1)]

    qd = quaddata(op, test_local_space, trial_local_space,
        [test_chart], refined_trial_chart, quadstrat)

    for (q,chart) in enumerate(refined_trial_chart)
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
        end end end
    end
end

function momintegrals_nested!(op::IntegralOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_chart, trial_chart, out, strat, quadstrat)

    momintegrals!(op, test_local_space, trial_local_space,
        test_chart, trial_chart, out, strat)
end

function momintegrals_trial_refines_test!(op::IntegralOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_chart, trial_chart, out, strat, quadstrat)

    momintegrals!(op, test_local_space, trial_local_space,
        test_chart, trial_chart, out, strat)
end


function momintegrals_trial_refines_test!(op::MWOperator3D,
    test_local_space::RTRefSpace, trial_local_space::RTRefSpace,
    test_chart, trial_chart, out, strat::SauterSchwabStrategy, quadstrat)

    # 1. Refine the test_chart
    p1, p2, p3 = test_chart.vertices

    # TODO: generalise this to include more general refinements
    e1 = cartesian(neighborhood(test_chart, (0,1/2)))
    e2 = cartesian(neighborhood(test_chart, (1/2,0)))
    e3 = cartesian(neighborhood(test_chart, (1/2,1/2)))

    ct = cartesian(center(test_chart))

    refined_test_chart = SA[
        simplex(ct, p1, e3),
        simplex(ct, e3, p2),
        simplex(ct, p2, e1),
        simplex(ct, e1, p3),
        simplex(ct, p3, e2),
        simplex(ct, e2, p1)]

    qd = quaddata(op, test_local_space, trial_local_space,
        refined_test_chart, [trial_chart], quadstrat)

    for (p,chart) in enumerate(refined_test_chart)
        qr = quadrule(op, test_local_space, trial_local_space,
            p, chart, 1, trial_chart, qd, quadstrat)

        Q = restrict(test_local_space, test_chart, chart)
        zlocal = zero(out)
        momintegrals!(op, test_local_space, trial_local_space,
            chart, trial_chart, zlocal, qr)

        for j in 1:3
            for i in 1:3
                for k in 1:3
                # out[i,j] += zlocal[i,k] * Q[j,k]
                out[i,j] += Q[i,k] * zlocal[k,j]
        end end end
    end
end