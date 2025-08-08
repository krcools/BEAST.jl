

abstract type Functional{T} end
scalartype(::Functional{T}) where {T} = T

defaultquadstrat(fn::Functional, basis) = SingleNumQStrat(8)
quaddata(fn::Functional, refs, cells, qs::SingleNumQStrat) =
    quadpoints(refs, cells, [qs.quad_rule])
quadrule(fn::Functional, refs, p, cell, qd, qs::SingleNumQStrat) = qd[1,p]

"""
    assemble(fn, tfs)

Assemble the vector of test coefficients corresponding to functional
`fn` and test functions `tfs`.
"""
function assemble(field::Functional, tfs;
    quadstrat=defaultquadstrat(field, tfs))

    R = scalartype(tfs)
    b = zeros(Complex{R}, numfunctions(tfs))
    store(v,m) = (b[m] += v)
    assemble!(field, tfs, store; quadstrat)
    return b
end

function assemble!(field::Functional, tfs::DirectProductSpace, store;
    quadstrat=defaultquadstrat(field, tfs))

    I = Int[0]
    for s in tfs.factors push!(I, last(I) + numfunctions(s)) end
    for (i,s) in enumerate(tfs.factors)
        store1(v,m) = store(v, m + I[i])
        assemble!(field, s, store1; quadstrat)
    end
end

function assemble!(field::Functional, tfs::Space, store;
    quadstrat=defaultquadstrat(field, tfs))

    qs = quadstrat(field, tfs)
    tels, tad = assemblydata(tfs)

    trefs = refspace(tfs)
    qd = quaddata(field, trefs, tels, qs)

    tgeo = geometry(tfs)
    tdom = domain(chart(tgeo, first(tgeo)))
    num_trefs = numfunctions(trefs, tdom)

    for (t, tcell) in enumerate(tels)

        # compute the testing with the reference elements
        qr = quadrule(field, trefs, t, tcell, qd, qs)
        blocal = celltestvalues(trefs, tcell, field, qr)

        for i in 1 : num_trefs
            for (m,a) in tad[t,i]
                store(a*blocal[i], m)
            end
        end

    end

end

function assemble!(field::Functional, tfs::subdBasis, store;
    quadstrat=defaultquadstrat(field, tfs))

    tels, tad = assemblydata(tfs)

    trefs = refspace(tfs)
    qd = quaddata(field, trefs, tels)

    for (t, tcell) in enumerate(tels)

        # compute the testing with the reference elements
        qr = quadrule(field, trefs, t, tcell, qd, quadstrat)
        blocal = celltestvalues(trefs, tcell, field, qr, quadstrat)

        for i in 1 : length(tad[t])
            for (m,a) in tad[t][i]
                store(a*blocal[i], m)
            end
        end

    end

end

function celltestvalues(tshs::RefSpace{T}, tcell, field, qr) where {T}

    num_tshs = numfunctions(tshs, domain(tcell))
    interactions = zeros(Complex{T}, num_tshs)

    num_oqp = length(qr)

    for p in 1 : num_oqp
        mp = qr[p].point

        dx =qr[p].weight

        fval = field(mp)
        tvals = qr[p].value

        for m in 1 : num_tshs
            tval = tvals[m]

            igd = integrand(field, tval, fval)
            interactions[m] += igd * dx
        end
    end

    return interactions
end

function celltestvalues(tshs::subReferenceSpace{T,D}, tcell, field, qr) where {T,D}

    num_oqp = length(qr)
    num_tshs = length(qr[1].value[1])
    interactions = (Complex{T}, num_tshs)
    for p in 1 : num_oqp
        mp = qr[p].point

        dx =qr[p].weight

        fval = field(mp)
        tvals = qr[p].value

        interactions = zeros(Complex{T}, num_tshs)
        for m in 1 : num_tshs
            tval = tvals[1][m]

            igd = integrand(field, tval, fval)
            interactions[m] += igd * dx
        end
    end

    return interactions
end



### Framework allowing to create Functional from Field by taking the trace
mutable struct CrossTraceMW{T,F} <: Functional{T}
    field::F
end

mutable struct TangTraceMW{T,F} <: Functional{T}
    field::F
end

CrossTraceMW(f::F) where {F} = CrossTraceMW{scalartype(f), F}(f)
CrossTraceMW{T}(f::F) where {T,F} = CrossTraceMW{T,F}(f)

TangTraceMW(f::F) where {F} = TangTraceMW{scalartype(f), F}(f)
TangTraceMW{T}(f::F) where {T,F} = TangTraceMW{T,F}(f)

# scalartype(x::CrossTraceMW) = scalartype(x.field)
# scalartype(::CrossTraceMW{T}) where {T} = T
# scalartype(::TangTraceMW{T}) where {T} = T


# scalartype(t::TangTraceMW) = scalartype(t.field)

cross(::NormalVector, p) = CrossTraceMW(p)
cross(t::CrossTraceMW, ::NormalVector) = TangTraceMW(t.field)

function (ϕ::CrossTraceMW)(p)
    F = ϕ.field
    x = cartesian(p)
    n = normal(p)
    return n × F(x)
end

function (ϕ::TangTraceMW)(p)
    F = ϕ.field
    x = cartesian(p)
    n = normal(p)
    return (n × F(x)) × n
end



struct NDotTrace{T,F} <: Functional{T}
    field::F
end


NDotTrace(f::F) where {F} = NDotTrace{scalartype(f), F}(f)
NDotTrace{T}(f::F) where {T,F} = NDotTrace{T,F}(f)
LinearAlgebra.dot(::NormalVector, f) = NDotTrace(f)

# scalartype(s::NDotTrace{T}) where {T} = T

(ϕ::NDotTrace)(p) = dot(normal(p), ϕ.field(cartesian(p)))

integrand(::Any, testvals, fieldval) = dot(testvals[1], fieldval)
# integrand(::TangTraceMW, gx, ϕx) = gx[1] ⋅ ϕx
# integrand(::CrossTraceMW, test_vals, field_val) = test_vals[1] ⋅ field_val
# integrand(::NDotTrace, g, ϕ) = dot(g.value, ϕ)