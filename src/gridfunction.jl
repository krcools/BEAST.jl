abstract type AbstractMeshFunction end

Base.:*(b::Number, fun::AbstractMeshFunction) = LinearCombinationOfAbstractMeshFunctions([b], [fun])

function Base.:+(fun1::AbstractMeshFunction, fun2::AbstractMeshFunction)
    if geometry(fun1) != geometry(fun2)
        error("Functions must be defined over the same geometry.")
    end
    return LinearCombinationOfAbstractMeshFunctions([1.0; 1.0], [fun1; fun2])
end

function Base.:-(fun1::AbstractMeshFunction, fun2::AbstractMeshFunction)
    if geometry(fun1) != geometry(fun2)
        error("Functions must be defined over the same geometry.")
    end
    return LinearCombinationOfAbstractMeshFunctions([1.0; -1.0], [fun1; fun2])
end

struct GlobalFunction{F,X,I} <: AbstractMeshFunction
    fun::F
    geo::X
    cells::Vector{I}
end

defaultquadstrat(field::GlobalFunction) = SingleNumQStrat(7)
quaddata(fn::GlobalFunction, simplexes, qs::SingleNumQStrat) = quadpoints(fn, simplexes, [qs.quad_rule])
quadrule(fn::GlobalFunction, refs, p, cell, qd, qs::SingleNumQStrat) = qd[1,p]
geometry(fn::GlobalFunction) = fn.geo

function (f::GlobalFunction{F, X, I})(mp::CompScienceMeshes.MeshPointNM) where {F, X, I}

    # Predict return type
	#PT = promote_type(coordtype(chart(mp)))

    return f.fun(cartesian(mp))
end

function (f::GlobalFunction{F, X, I})(x) where {F, X, I}
    return f.fun(x)
end

struct FEMFunction{T, X} <: AbstractMeshFunction
    coeffs::AbstractVector{T}
    space::X
end

function Base.getindex(u::BEAST.FEMFunction, b::BEAST.Variational.HilbertVector)
    return BEAST.FEMFunction(u.coeffs[b], u.space[b.idx])
end

function Base.:-(u::BEAST.FEMFunction)
    return BEAST.FEMFunction(-u.coeffs, u.space)
end

@testitem "FEMFunction getindex" begin
    using CompScienceMeshes

    Γ = CompScienceMeshes.meshrectangle(1.0, 1.0, 0.5)
    X1 = BEAST.raviartthomas(Γ)

    X = X1 × X1
    ax = BEAST.NestedUnitRanges.nestedrange(X, 1, numfunctions)

    coeffs = rand(numfunctions(X))
    coeffs = BEAST.BlockArrays.BlockedVector(coeffs, (ax,))

    u = BEAST.FEMFunction(coeffs, X)

    @hilbertspace m[1:2]
    u1 = u[m[1]]
    u2 = u[m[2]]
    u2 = -u1

    @test u1 isa BEAST.FEMFunction
    @test u2 isa BEAST.FEMFunction

    @test u1.coeffs == -u2.coeffs
end

defaultquadstrat(field::FEMFunction) = SingleNumQStrat(7)
returntype(f::FEMFunction{T,X}) where {T,X} = T

quaddata(fn::FEMFunction, refs, cells, qs::SingleNumQStrat) = quadpoints(refs, cells, [qs.quad_rule])
quadrule(fn::FEMFunction, refs, p, cell, qd, qs::SingleNumQStrat) = qd[1,p]
geometry(fn::FEMFunction) = geometry(fn.space)

mutable struct LinearCombinationOfAbstractMeshFunctions{T}
    coeffs::Vector{T}
    gfs::Vector
end

geometry(lincombfun::LinearCombinationOfAbstractMeshFunctions) = geometry(first(lincombfun.gfs))

Base.length(lc::LinearCombinationOfAbstractMeshFunctions) = length(lc.coeffs)

Base.:*(b::Number, lcb::LinearCombinationOfAbstractMeshFunctions) =
    LinearCombinationOfAbstractMeshFunctions(b .* lcb.coeffs, lcb.gfs)

function Base.:+(lcb::LinearCombinationOfAbstractMeshFunctions{T}, fun::AbstractMeshFunction) where {T}

    if geometry(lcb) != geometry(fun)
        error("Functions must be defined over the same geometry.")
    end

    return LinearCombinationOfAbstractMeshFunctions([lcb.coeffs[:]; T(1)], [lcb.gfs[:]; fun])
end

function Base.:-(lcb::LinearCombinationOfAbstractMeshFunctions{T}, fun::AbstractMeshFunction) where {T}

    if geometry(lcb) != geometry(fun)
        error("Functions must be defined over the same geometry.")
    end

    return LinearCombinationOfAbstractMeshFunctions([lcb.coeffs[:]; -T(1)], [lcb.gfs[:]; fun])
end

Base.:+(fun::AbstractMeshFunction, lcb::LinearCombinationOfAbstractMeshFunctions) = lcb + fun

Base.:-(fun::AbstractMeshFunction, lcb::LinearCombinationOfAbstractMeshFunctions) = lcb - fun

function Base.:+(lcb1::LinearCombinationOfAbstractMeshFunctions, lcb2::LinearCombinationOfAbstractMeshFunctions)

    if geometry(lcb1) != geometry(lcb2)
        error("Functions must be defined over the same geometry.")
    end

    return LinearCombinationOfAbstractMeshFunctions([lcb1.coeffs[:]; lcb2.coeffs[:]], [lcb1.gfs[:]; lcb2.gfs[:]])
end

function Base.:-(lcb1::LinearCombinationOfAbstractMeshFunctions{T1}, lcb2::LinearCombinationOfAbstractMeshFunctions{T2}) where {T1, T2}

    if geometry(lcb1) != geometry(lcb2)
        error("Functions must be defined over the same geometry.")
    end

    return LinearCombinationOfAbstractMeshFunctions([lcb1.coeffs[:]; -T2(1) .* lcb2.coeffs[:]], [lcb1.gfs[:]; lcb2.gfs[:]])
end

Base.iterate(lc::LinearCombinationOfAbstractMeshFunctions, state=1) =
    state > length(lc) ? nothing : ((lc.coeffs[state], lc.gfs[state]), state+1)

defaultquadstrat(field::LinearCombinationOfAbstractMeshFunctions) = SingleNumQStrat(7)

"""
    Lp_integrate

Integrates FEMFunctions and GlobalFunctions as well as linear combinations.
NOTE: FEMFunctions on the dual mesh are currently not supported.
"""

function Lp_integrate(gf::FEMFunction; p=2, quadstrat=defaultquadstrat(gf))
    return Lp_integrate(LinearCombinationOfAbstractMeshFunctions([1.0], [gf]); p=p, quadstrat=quadstrat)
end

function Lp_integrate(
    glf::GlobalFunction{F, X, I};
    p=2,
    quadstrat=defaultquadstrat(glf)
) where {F, X, I}
    return Lp_integrate(LinearCombinationOfAbstractMeshFunctions([1.0], [glf]); p=p, quadstrat=quadstrat)
end

"""
    indices_splitfemglobal(lincombgfs)

Given a linear combination of FEMFunctions and/or GlobalFunctions, return the indices
of the fem and the global functions
"""
function indices_splitfemglobal(lincombgfs::LinearCombinationOfAbstractMeshFunctions{T}) where T

    indices_fem = Int64[]
    indices_global = Int64[]

    for (i, gf) in enumerate(lincombgfs.gfs)
        if gf isa FEMFunction
            push!(indices_fem, i)
        else
            push!(indices_global, i)
        end
    end

    return indices_fem, indices_global
end

function Lp_integrate(
    lincombgfs::LinearCombinationOfAbstractMeshFunctions{T};
    p=2,
    quadstrat=SingleNumQStrat(7)
) where {T}

    # TODO: Add routine that geos are consistent
    indices_fem, indices_global = indices_splitfemglobal(lincombgfs)

    if length(indices_fem) >= 1
        geo = geometry(lincombgfs.gfs[indices_fem[1]].space)
    elseif length(indices_global) >= 1
        geo = lincombgfs.gfs[indices_global[1]].geo
    else
        error("Error: empty linear combination")
    end

    asmdata = [ assemblydata(lincombgfs.gfs[i].space, onlyactives=false) for i in indices_fem]

    qds_fem = Any[]
    trefs_fem = Any[]

    numquadpoints = 0
    for (i, j) in enumerate(indices_fem)
        push!(trefs_fem, refspace(lincombgfs.gfs[j].space))
        tels = asmdata[i][1]
        qd = quaddata(lincombgfs.gfs[j], trefs_fem[i], tels, quadstrat)
        push!(qds_fem, qd)
    end

    qds_global = Any[]
    active_global = Any[]
    for (i, j) in enumerate(indices_global)
        push!(active_global, falses(numcells(lincombgfs.gfs[j].geo)))

        for cellindex in lincombgfs.gfs[j].cells
            active_global[i][cellindex] = true
        end

        tels = instantiate_charts(lincombgfs.gfs[j].geo, numcells(lincombgfs.gfs[j].geo), trues(numcells(lincombgfs.gfs[j].geo)))
        qd = quaddata(lincombgfs.gfs[j], tels, quadstrat)
        push!(qds_global, qd)
    end

    if length(indices_fem) >= 1
        numquadpoints = length(qds_fem[1][1, 1])
    elseif length(indices_global) >= 1
        numquadpoints = length(qds_global[1][1, 1])
    else
        error("Error: empty linear combination")
    end

    # It is difficult to predict the return type of an arbitrary function
    quadrature_sampling_gf_res = zeros(ComplexF64, numquadpoints, numcells(geo))

    retval = Float64(0) #returntype(first(lincombgfs.gfs))(0)

    tels = instantiate_charts(geo, numcells(geo), trues(numcells(geo)))

    for (t, tcell) in enumerate(tels)

        qrs_fem = [quadrule(lincombgfs.gfs[j], trefs_fem[i], t, tcell, qds_fem[i], quadstrat) for (i, j) in enumerate(indices_fem)]

        # Loop over FEM functions first

        #for (i, (~, gftad, ~)) in enumerate(asmdata)
        for (i, j) in enumerate(indices_fem)
            ~, gftad, ~ = asmdata[i]
            cellquadsamplingvaluesvalues!(
                quadrature_sampling_gf_res,
                t,
                trefs_fem[i],
                gftad,
                lincombgfs.coeffs[j],
                lincombgfs.gfs[j],
                qrs_fem[i]
            )
        end

        qrs_global = [
            quadrule(lincombgfs.gfs[j], nothing, t, tcell, qds_global[i], quadstrat)
            for (i, j) in enumerate(indices_global)
        ]

        # TODO: loop over grid functions

        for (i, j) in enumerate(indices_global)

            # We only integrate, if the cell is part of the domain
            active_global[i][t] == false && continue

            cellquadsamplingvaluesvalues!(
                quadrature_sampling_gf_res,
                t,
                nothing,
                nothing,
                lincombgfs.coeffs[j],
                lincombgfs.gfs[j],
                qrs_global[i]
            )
        end

        # The weights are all the same since we use the same quadrature
        # Just pick any from the no-empty sets
        if length(indices_fem) >= 1
            qr = first(qrs_fem)
        elseif length(indices_global) >= 1
            qr = first(qrs_global)
        else
            error("Error: empty linear combination")
        end

        for i = 1:numquadpoints
            retval += abs(quadrature_sampling_gf_res[i, t])^p*qr[i].weight
        end
    end

    return retval^(1/p)
end

function cellquadsamplingvaluesvalues!(
    quadsampling,
    cellindex,
    tshs::RefSpace,
    tad,
    lincombcoeff,
    gf,
    qr
)

    # num_tshs = numfunctions(tshs)

    num_oqp = length(qr)
    
    data = tad.data

    for p in 1 : num_oqp

        tvals = qr[p].value

        for m in 1 : length(tvals)
        # for m in 1 : num_tshs
            tval = tvals[m]

            for k = 1:size(data, 1)
                # w is the is weight of shape function
                # of the global basis function, omitted via ~
                b, w = data[k, m, cellindex] 

                igd = tval[1]*w*lincombcoeff*gf.coeffs[b]
                quadsampling[p, cellindex] += igd 
            end
        end
    end

end

function cellquadsamplingvaluesvalues!(
    quadsampling,
    cellindex,
    tshs,
    tad::Nothing,
    lincombcoeff,
    gf::GlobalFunction,
    qr
)

    num_oqp = length(qr)

    for p in 1 : num_oqp
        igd = qr[p].value
        if lincombcoeff isa Number
            igd *= lincombcoeff
        end
        quadsampling[p, cellindex] += igd 
    end

end