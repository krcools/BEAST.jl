struct Integrand{Op,LSt,LSb,Elt,Elb}
    operator::Op
    local_test_space::LSt
    local_trial_space::LSb
    test_chart::Elt
    trial_chart::Elb
end


function (igd::Integrand)(u,v)
    
    x = neighborhood(igd.test_chart,u)
    y = neighborhood(igd.trial_chart,v)
    
    f = igd.local_test_space(x)
    g = igd.local_trial_space(y)

    return jacobian(x) * jacobian(y) * igd(x,y,f,g)
end

# For divergence conforming basis and trial functions, an alternative evaluation
# of the integrand is possible that avoids the computation of the chart jacobian
# determinants.
function (igd::Integrand{<:IntegralOperator,<:DivRefSpace,<:DivRefSpace})(u,v)
    test_domain = CompScienceMeshes.domain(igd.test_chart)
    bsis_domain = CompScienceMeshes.domain(igd.trial_chart)

    x = CompScienceMeshes.neighborhood_lazy(igd.test_chart,u)
    y = CompScienceMeshes.neighborhood_lazy(igd.trial_chart,v)
    
    p = neighborhood(test_domain, u)
    q = neighborhood(bsis_domain, v)
    
    f̂ = igd.local_test_space(p)
    ĝ = igd.local_trial_space(q)

    Dx = tangents(x)
    Dy = tangents(y)

    f = map(f̂) do fi
        (value = Dx * fi.value, divergence = fi.divergence) end
    g = map(ĝ) do gi
        (value = Dy * gi.value, divergence = gi.divergence) end

    igd(x,y,f,g)
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




function _integrands_leg_gen(f::Type{U}, g::Type{V}) where {U<:SVector{N}, V<:SVector{M}} where {M,N}
    ex = :(SMatrix{N,M}(()))
    for m in 1:M
        for n in 1:N
            push!(ex.args[2].args, :(integrand(op, kervals, f[$n], x, g[$m], y)))
        end
    end
    return ex
end
@generated function _integrands_leg(op, kervals, f::SVector{N}, x, g::SVector{M}, y) where {M,N}
    _integrands_leg_gen(f, g)
end


# Support for legacy kernels
function (igd::Integrand)(x,y,f,g)

    op = igd.operator
    kervals = kernelvals(op, x, y)
    _integrands_leg(op, kervals, f, x, g, y)

end

# function CompScienceMeshes.permute_vertices(
#     ch::CompScienceMeshes.RefQuadrilateral, I)

#     V = vertices(ch)
#     return Quadrilateral(V[I[1]], V[I[2]], V[I[3]], V[I[4]])
# end

struct PulledBackIntegrand{I,C1,C2}
    igd::I
    chart1::C1
    chart2::C2
end

function (f::PulledBackIntegrand)(u,v)
    # In general I think a Jacobian determinant needs to be included. For Simplical and
    # Quadrilateral charts this is not needed because they are 1.
    f.igd(cartesian(f.chart1,u), cartesian(f.chart2,v))
end

function pulledback_integrand(igd,
    I, chart1,
    J, chart2)

    dom1 = domain(chart1)
    dom2 = domain(chart2)

    ichart1 = CompScienceMeshes.permute_vertices(dom1, I)
    ichart2 = CompScienceMeshes.permute_vertices(dom2, J)

    PulledBackIntegrand(igd, ichart1, ichart2)
end 

function momintegrals!(op::Operator,
    test_local_space, trial_local_space,
    test_chart, trial_chart,
    out, rule::SauterSchwabStrategy)

    I, J, _, _ = SauterSchwabQuadrature.reorder(
        vertices(test_chart),
        vertices(trial_chart), rule)

    num_tshapes = numfunctions(test_local_space, domain(test_chart))
    num_bshapes = numfunctions(trial_local_space, domain(trial_chart))

    igd = Integrand(op, test_local_space, trial_local_space, test_chart, trial_chart)
    igdp = pulledback_integrand(igd, I, test_chart, J, trial_chart)
    G = SauterSchwabQuadrature.sauterschwab_parameterized(igdp, rule)
    out[1:num_tshapes, 1:num_bshapes] .+= G

    nothing
end



