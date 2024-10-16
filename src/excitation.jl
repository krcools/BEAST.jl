

abstract type Functional end

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

    tels, tad = assemblydata(tfs)

    trefs = refspace(tfs)
    qd = quaddata(field, trefs, tels, quadstrat)

    tgeo = geometry(tfs)
    tdom = domain(chart(tgeo, first(tgeo)))
    num_trefs = numfunctions(trefs, tdom)

    for (t, tcell) in enumerate(tels)

        # compute the testing with the reference elements
        qr = quadrule(field, trefs, t, tcell, qd, quadstrat)
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
