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

function sauterschwab_parameterized(igdp, rule::SauterSchwabStrategy)
    return SauterSchwabQuadrature.sauterschwab_parameterized(igdp, rule)
end

function sauterschwab_parameterized(igdp, rule::SauterSchwabQuadrature1D.SauterSchwabStrategy1D)
    return SauterSchwabQuadrature1D.sauterschwab_parameterized1D(igdp, rule)
end

function sauterschwab_reorder(test_vertices, trial_vertices, rule::SauterSchwabStrategy)
    I, J, _, _ = SauterSchwabQuadrature.reorder(test_vertices, trial_vertices, rule)

    return I, J
end

function sauterschwab_reorder(test_vertices, trial_vertices, rule::SauterSchwabQuadrature1D.SauterSchwabStrategy1D)
    I, J, _, _ = SauterSchwabQuadrature1D.reorder(test_vertices, trial_vertices, rule)

    return I, J
end

function momintegrals!(op::Operator,
    test_local_space, trial_local_space,
    test_chart, trial_chart,
    out, rule::Union{SauterSchwabStrategy,SauterSchwabQuadrature1D.SauterSchwabStrategy1D})

    I, J = sauterschwab_reorder(
        vertices(test_chart),
        vertices(trial_chart),
        rule
    )

    num_tshapes = numfunctions(test_local_space, domain(test_chart))
    num_bshapes = numfunctions(trial_local_space, domain(trial_chart))

    igd = Integrand(op, test_local_space, trial_local_space, test_chart, trial_chart)
    igdp = pulledback_integrand(igd, I, test_chart, J, trial_chart)

    G = sauterschwab_parameterized(igdp, rule)

    for j in 1:num_bshapes
        for i in 1:num_tshapes
            out[i,j] += G[i,j]
        end
    end

    nothing
end

function momintegrals!(op::Operator,
    test_local_space, trial_local_space,
    test_chart, trial_chart,
    out, rule::SauterSchwab3DStrategy)

    I, J = SauterSchwab3D.reorder(rule.sing)

    num_tshapes = numfunctions(test_local_space, domain(test_chart))
    num_bshapes = numfunctions(trial_local_space, domain(trial_chart))

    igd = Integrand(op, test_local_space, trial_local_space, test_chart, trial_chart)
    igdp = pulledback_integrand(igd, I, test_chart, J, trial_chart)
    G = SauterSchwab3D.sauterschwab_parameterized(igdp, rule)
    out[1:num_tshapes, 1:num_bshapes] .+= G
    nothing
end

# abstract type Singularity end
# abstract type Singularity6D <: Singularity end
# abstract type Singularity5D <: Singularity end
# abstract type Singularity4D <: Singularity end

# struct Singularity6DPositiveDistance  <: Singularity6D end
# struct Singularity6DPoint  <: Singularity6D T::SVector{1,Int64}; S::SVector{1,Int64} end
# struct Singularity6DEdge   <: Singularity6D T::SVector{2,Int64}; S::SVector{2,Int64} end
# struct Singularity6DFace   <: Singularity6D T::SVector{3,Int64}; S::SVector{3,Int64} end
# struct Singularity6DVolume <: Singularity6D T::SVector{4,Int64}; S::SVector{4,Int64} end

# struct Singularity5DPositiveDistance  <: Singularity5D end
# struct Singularity5DPoint  <: Singularity5D T::SVector{1,Int64}; S::SVector{1,Int64} end
# struct Singularity5DEdge   <: Singularity5D T::SVector{2,Int64}; S::SVector{2,Int64} end
# struct Singularity5DFace   <: Singularity5D T::SVector{3,Int64}; S::SVector{3,Int64} end

# struct Singularity4DPositiveDistance  <: Singularity4D end
# struct Singularity4DPoint  <: Singularity4D T::SVector{1,Int64}; S::SVector{1,Int64} end
# struct Singularity4DEdge   <: Singularity4D T::SVector{2,Int64}; S::SVector{2,Int64} end
# struct Singularity4DFace   <: Singularity4D T::SVector{3,Int64}; S::SVector{3,Int64} end

reversestrat(a::SauterSchwab3D.Singularity6DPoint) = SauterSchwab3D.Singularity6DPoint(a.S, a.T)
reversestrat(a::SauterSchwab3D.Singularity6DEdge) = SauterSchwab3D.Singularity6DEdge(a.S, a.T)
reversestrat(a::SauterSchwab3D.Singularity6DFace) = SauterSchwab3D.Singularity6DFace(a.S, a.T)
reversestrat(a::SauterSchwab3D.Singularity6DVolume) = SauterSchwab3D.Singularity6DVolume(a.S, a.T)

reversestrat(a::SauterSchwab3D.Singularity5DPoint) = SauterSchwab3D.Singularity5DPoint(a.S, a.T)
reversestrat(a::SauterSchwab3D.Singularity5DEdge) = SauterSchwab3D.Singularity5DEdge(a.S, a.T)
reversestrat(a::SauterSchwab3D.Singularity5DFace) = SauterSchwab3D.Singularity5DFace(a.S, a.T)

reversestrat(a::SauterSchwab3D.Singularity4DPoint) = SauterSchwab3D.Singularity4DPoint(a.S, a.T)
reversestrat(a::SauterSchwab3D.Singularity4DEdge) = SauterSchwab3D.Singularity4DEdge(a.S, a.T)
reversestrat(a::SauterSchwab3D.Singularity4DFace) = SauterSchwab3D.Singularity4DFace(a.S, a.T)
reversestrat(a::T) where {T <: SauterSchwab3D.SauterSchwab3DStrategy} = T(reversestrat(a.sing),a.qps)

function momintegrals!(op::Operator,
    test_local_space, trial_local_space,
    test_chart, trial_chart,
    out, rule::_TransposedStrat{<:SauterSchwab3DStrategy})
    rule2 = reversestrat(rule.strat)
    J, I = SauterSchwab3D.reorder(rule2.sing)
    #I2,J2 = SauterSchwab3D.reorder(rule.strat.sing)
    #@assert I == I2
    num_tshapes = numfunctions(test_local_space, domain(test_chart))
    num_bshapes = numfunctions(trial_local_space, domain(trial_chart))

    igd = (u,v) -> Integrand(op, test_local_space, trial_local_space, test_chart, trial_chart)(v,u)
    igdp = pulledback_integrand(igd, J, trial_chart, I, test_chart)
    G = SauterSchwab3D.sauterschwab_parameterized(igdp, rule2)
    out[1:num_tshapes, 1:num_bshapes] .+= G
    nothing
end