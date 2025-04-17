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
    return _integrands_leg(op, kervals, f, x, g, y)

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
    # out[1:num_tshapes, 1:num_bshapes] .+= G
    for j in 1:num_bshapes
        for i in 1:num_tshapes
            out[i,j] += G[i,j]
        end
    end

    nothing
end

struct Integrand1D{Op,LSt,LSb,Elt,Elb}
    operator::Op
    local_test_space::LSt
    local_trial_space::LSb
    test_chart::Elt
    trial_chart::Elb
end

function (igd::Integrand1D)(x,y,f,g)

    op = igd.operator

    kervals = kernelvals(op, x, y)

    wynik=_integrands_leg(op, kervals, f, x, g, y)
    return wynik
end

function (igd::Integrand1D)(u,v)

    x = neighborhood(igd.test_chart,u)
    y = neighborhood(igd.trial_chart,v)

    f = igd.local_test_space(x)
    g = igd.local_trial_space(y)

    return jacobian(x) * jacobian(y) * igd(x,y,f,g)
end

struct PulledBackIntegrand1D{I,C}
    igd::I
    chart1::C
    chart2::C
end

function (f::PulledBackIntegrand1D)(u,v)

    return f.igd(u, v)
end

function reference_vertices(::Type{CompScienceMeshes.ReferenceSimplex{1, T, 2}}) where T
    return (SVector{2,T}(0,0), SVector{2,T}(1,0))
end

#=
function pulledback_integrand1D(igd,
    I, chart1,
    J, chart2)

    dom1 = domain(chart1)
    dom2 = domain(chart2)

    V1 = reference_vertices(typeof(dom1))
    V2 = reference_vertices(typeof(dom2))

    ichart1 = simplex(V1[I[1]], V1[I[2]])
    ichart2 = simplex(V2[J[1]], V2[J[2]])

    PulledBackIntegrand1D(igd, ichart1, ichart2)
end
=#

function pulledback_integrand1D(igd,
    chart1,
    chart2)
    PulledBackIntegrand1D(igd, chart1, chart2)
end

function momintegrals!(op::Operator,
    test_local_space, trial_local_space,
    test_chart, trial_chart,
    out, rule::BEAST.SauterSchwabQuadrature1D.CommonEdge)


    igd = Integrand1D(op, test_local_space, trial_local_space, test_chart, trial_chart)

    G = BEAST.SauterSchwabQuadrature1D.sauterschwab_parameterized1D(igd, rule)
    out[
        1:numfunctions(test_local_space, domain(test_chart)),
        1:numfunctions(trial_local_space, domain(trial_chart))
    ] .+= G
    nothing
end


"""
it is the test for correctly order of vertices, which results in anticlock way r1->r3->r2->r1
this means that for u,v = 0 r1, r3 and for u,v=1 the expected common node r2 is present
if this condition is not met, a reordering function must be implemented
"""

function reverse_chart(chart::CompScienceMeshes.Simplex{2,1,1,2,T}) where T
    new_verts = SVector{2,SVector{2,T}}(chart.vertices[2], chart.vertices[1])

    new_tangent = -chart.tangents[1]

    # TODO: double check normal should not flip
    new_normal = chart.normals[1]

    new_vol = chart.volume

    return CompScienceMeshes.Simplex{2,1,1,2,T}(
        new_verts,
        SVector{1,SVector{2,T}}(new_tangent),
        SVector{1,SVector{2,T}}(new_normal),
        new_vol
    )
end


function momintegrals!(op::Operator,
    test_local_space, trial_local_space,
    test_chart, trial_chart,
    out, rule::BEAST.SauterSchwabQuadrature1D.CommonVertex)

    igd = Integrand1D(op, test_local_space, trial_local_space, test_chart, trial_chart)

    u = 0.0
    v = 1.0

    if !(cartesian(igd.test_chart, u) ≈ cartesian(igd.trial_chart, v))

        itest_chart = reverse_chart(igd.test_chart)
        itrial_chart = reverse_chart(igd.trial_chart)

        igdp = pulledback_integrand1D(igd, itest_chart, itrial_chart)

        G = BEAST.SauterSchwabQuadrature1D.sauterschwab_parameterized1D(igdp, rule)
        out[
            1:numfunctions(test_local_space, domain(test_chart)),
            1:numfunctions(trial_local_space, domain(trial_chart))
        ] .+= G

    else

        G = BEAST.SauterSchwabQuadrature1D.sauterschwab_parameterized1D(igd, rule)
        out[
            1:numfunctions(test_local_space, domain(test_chart)),
            1:numfunctions(trial_local_space, domain(trial_chart))
        ] .+= G
    end
    nothing
end

