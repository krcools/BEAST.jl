export assemble

abstract type Functional end

quaddata(fn::Functional, refs, cells) = quadpoints(refs, cells, [5])
quadrule(fn::Functional, refs, p, cell, qd) = qd[1,p]

"""
    assemble(fn, tfs)

Assemble the vector of test coefficients corresponding to functional
`fn` and test functions `tfs`.
"""
function assemble(field::Functional, tfs)

    b = zeros(Complex128, numfunctions(tfs))
    store(v,m) = (b[m] += v)
    assemble!(field, tfs, store)
    return b
end

function assemble!(field::Functional, tfs::DirectProductSpace, store)
    I = Int[0]
    for s in tfs.factors push!(I, last(I) + numfunctions(s)) end
    for (i,s) in enumerate(tfs.factors)
        store1(v,m) = store(v, m + I[i])
        assemble!(field, s, store1)
    end
end

function assemble!(field::Functional, tfs::Space, store)

    tels, tad = assemblydata(tfs)

    trefs = refspace(tfs)
    qd = quaddata(field, trefs, tels)

    for (t, tcell) in enumerate(tels)

        # compute the testing with the reference elements
        qr = quadrule(field, trefs, t, tcell, qd)
        blocal = celltestvalues(trefs, tcell, field, qr)

        for i in 1 : numfunctions(trefs)
            for (m,a) in tad[t,i]
                store(a*blocal[i], m)
            end
        end

    end

end

function celltestvalues(tshs, tcell, field, qr)

    num_tshs = numfunctions(tshs)
    interactions = zeros(Complex128, num_tshs)

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
