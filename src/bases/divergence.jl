export divergence

function divergence(X::DirectProductSpace{T}) where T
    x = Space{T}[divergence(s) for s in X.factors]
    DirectProductSpace(x)
end

"""
    divergence(x)

Compute the divergence of a finite element space.
"""
function divergence(x::Space)
    ref = refspace(x)
    geo = geometry(x)
    els = elements(geo)
    fns = x.fns
    dvs = similar(fns)
    for (i,fn) in enumerate(fns)
        dvs[i] = similar(fns[i])
        for (j,sh) in enumerate(fn)
            el = els[sh.cellid]
            dvs[i][j] = divergence(ref, sh, el)
        end
    end
    divergence(x, geo, dvs)
end


"""
    curl(X)

Compute the curl of a finite element basis. The resulting set of functions might be linearly dependent because of the kernel of the curl operator.
"""
function curl(x::Space)
    ref = refspace(x)
    geo = geometry(x)
    els = elements(geo)
    crl = similar(x.fns)
    for (i,fn) in enumerate(x.fns)
        crl[i] = similar(x.fns[i], 0)
        for (j,sh) in enumerate(fn)
            el = els[sh.cellid]
            append!(crl[i], curl(ref, sh, el))
        end
    end
    curl(x, geo, crl)
end
