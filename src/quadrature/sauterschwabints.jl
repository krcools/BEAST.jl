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

# TODO: Remove this when the next version of CSM is relased!
function CompScienceMeshes.permute_vertices(
    ch::CompScienceMeshes.RefQuadrilateral, I)

    V = vertices(ch)
    return Quadrilateral(V[I[1]], V[I[2]], V[I[3]], V[I[4]])
end

struct PulledBackIntegrand{I,C}
    igd::I
    chart1::C
    chart2::C
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

    igd = Integrand(op, test_local_space, trial_local_space, test_chart, trial_chart)
    igdp = pulledback_integrand(igd, I, test_chart, J, trial_chart)
    G = SauterSchwabQuadrature.sauterschwab_parameterized(igdp, rule)
    out[1:numfunctions(test_local_space),1:numfunctions(trial_local_space)] .+= G

    nothing
end


# function momintegrals_test_refines_trial!(out, op,
#     test_functions, test_cell, test_chart,
#     trial_functions, trial_cell, trial_chart,
#     quadrule, quadstrat)

#     # test_local_space = refspace(test_functions)
#     # trial_local_space = refspace(trial_functions)

#     momintegrals!(out, op,
#         test_functions, test_cell, test_chart,
#         trial_functions, trial_cell, trial_chart, quadrule)
# end

# # const MWOperator3D = Union{MWSingleLayer3D, MWDoubleLayer3D}
# function momintegrals_test_refines_trial!(out, op,
#     test_functions, test_cell, test_chart,
#     trial_functions, trial_cell, trial_chart,
#     qr::SauterSchwabStrategy, quadstrat)

#     test_local_space = refspace(test_functions)
#     trial_local_space = refspace(trial_functions)

#     test_mesh = geometry(test_functions)
#     # trial_mesh = geometry(trial_functions)

#     parent_mesh = CompScienceMeshes.parent(test_mesh)
#     trial_charts = [chart(test_mesh, p) for p in CompScienceMeshes.children(parent_mesh, trial_cell)]

#     qd = quaddata(op, test_local_space, trial_local_space,
#         [test_chart], trial_charts, quadstrat)

#     for (q,chart) in enumerate(trial_charts)
#         qr = quadrule(op, test_local_space, trial_local_space,
#             1, test_chart, q ,chart, qd, quadstrat)
#         # @show qr

#         Q = restrict(trial_local_space, trial_chart, chart)
#         zlocal = zero(out)

#         momintegrals!(zlocal, op,
#             test_functions, test_cell, test_chart,
#             trial_functions, nothing, chart, qr)

#         for j in 1:numfunctions(trial_local_space)
#             for i in 1:numfunctions(test_local_space)
#                 for k in 1:size(Q, 2)
#                     out[i,j] += zlocal[i,k] * Q[j,k]
# end end end end end



# function momintegrals_trial_refines_test!(out, op,
#     test_functions, test_cell, test_chart,
#     trial_functions, trial_cell, trial_chart,
#     quadrule, quadstrat)

#     test_local_space = refspace(test_functions)
#     trial_local_space = refspace(trial_functions)

#     momintegrals!(out, op,
#         test_functions, test_cell, test_chart,
#         trial_functions, trial_cell, trial_chart, quadrule)
# end


# function momintegrals_trial_refines_test!(out, op,
#     test_functions, test_cell, test_chart,
#     trial_functions, trial_cell, trial_chart,
#     qr::SauterSchwabStrategy, quadstrat)

#     test_local_space = refspace(test_functions)
#     trial_local_space = refspace(trial_functions)

#     test_mesh = geometry(test_functions)
#     trial_mesh = geometry(trial_functions)

#     parent_mesh = CompScienceMeshes.parent(trial_mesh)
#     test_charts = [chart(trial_mesh, p) for p in CompScienceMeshes.children(parent_mesh, test_cell)] 

#     qd = quaddata(op, test_local_space, trial_local_space,
#         test_charts, [trial_chart], quadstrat)

#     for (p,chart) in enumerate(test_charts)
#         qr = quadrule(op, test_local_space, trial_local_space,
#             p, chart, 1, trial_chart, qd, quadstrat)

#         Q = restrict(test_local_space, test_chart, chart)
#         zlocal = zero(out)
#         momintegrals!(zlocal, op,
#             test_functions, nothing, chart,
#             trial_functions, trial_cell, trial_chart, qr)

#         for j in 1:numfunctions(trial_local_space)
#             for i in 1:numfunctions(test_local_space)
#                 for k in 1:size(Q, 2)
#                     out[i,j] += Q[i,k] * zlocal[k,j]
#         end end end
#     end
# end
