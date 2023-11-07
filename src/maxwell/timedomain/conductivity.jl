abstract type ConductivityTD_Functionaltype <: BEAST.Functional end

abstract type ConductivityTD_Operatortype <: BEAST.LocalOperator end

mutable struct ConductivityTDFunc <: ConductivityTD_Functionaltype
    chr::BEAST.Polynomial
    dchr::BEAST.Polynomial
    numdiffs::Int
    efield::Array{SVector{3, Float64},2}
    jflux::Array{SVector{3, Float64},2}
end

mutable struct ConductivityTDOp <: ConductivityTD_Operatortype
    op::ConductivityTDFunc
end

scalartype(p::ConductivityTDFunc) = eltype(p.efield[1])
scalartype(p::ConductivityTDOp) = eltype(p.op.efield[1])

function conductivityfunc(chr::BEAST.Polynomial;numdiffs=0)
    dchr = derive(chr)
    efield = fill(SVector(0.0,0.0,0.0), (2,2))
    ConductivityTDFunc(chr, dchr, numdiffs, efield, efield)
end

function (f::ConductivityTDFunc)(cell, cqdpt, mp)
    ei = f.efield[cell, cqdpt]
    return f.chr(norm(ei))*ei
end

BEAST.integrand(::BEAST.ConductivityTD_Functionaltype, gx, ϕx) = gx[1] ⋅ ϕx
BEAST.defaultquadstrat(::BEAST.ConductivityTD_Functionaltype, ::BEAST.LinearRefSpaceTriangle, ::BEAST.LinearRefSpaceTriangle) = BEAST.SingleNumQStrat(10)
BEAST.defaultquadstrat(::BEAST.ConductivityTD_Functionaltype, ::BEAST.Space) = BEAST.SingleNumQStrat(10)

function quaddata(exc::BEAST.ConductivityTD_Functionaltype, g::BEAST.LinearRefSpaceTriangle, f::BEAST.LinearRefSpaceTriangle, tels, bels,
    qs::BEAST.SingleNumQStrat)

    u, w = BEAST.trgauss(qs.quad_rule)
    qd = [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
    A = BEAST._alloc_workspace(qd, g, f, tels, bels)

    return qd, A
end

function quadrule(exc::BEAST.ConductivityTD_Functionaltype, ψ::BEAST.RefSpace, ϕ::BEAST.RefSpace, τ, (qd,A), qs::BEAST.SingleNumQStrat)
    for i in eachindex(qd)
        q = qd[i]
        w, p = q[1], BEAST.neighborhood(τ,q[2])
        A[i] = (w, p, ψ(p), ϕ(p))
    end
    return A
end

function kernelvals(f::ConductivityTDOp, mp, cell, cqdpt)
    ei = f.op.efield[cell, cqdpt]
    if norm(ei)==0
        ei = 1e-9.+ei
    end
    dsigma = f.op.dchr(norm(ei))*kron(ei, ei')/norm(ei)+(f.op.chr(norm(ei)))*I(3)
    return dsigma
end

BEAST.integrand(op::ConductivityTD_Operatortype, kernel, x, g, f) = dot(g[1],kernel*f[1])

function assemble_local_matched!(biop::ConductivityTDOp, tfs::BEAST.Space, bfs::BEAST.Space, store;
    quadstrat=BEAST.defaultquadstrat(biop, tfs, bfs))

    tels, tad, ta2g = BEAST.assemblydata(tfs)
    bels, bad, ba2g = BEAST.assemblydata(bfs)

    bg2a = zeros(Int, length(BEAST.geometry(bfs)))
    for (i,j) in enumerate(ba2g) bg2a[j] = i end

    trefs = BEAST.refspace(tfs)
    brefs = BEAST.refspace(bfs)

    qd = BEAST.quaddata(biop, trefs, brefs, tels, bels, quadstrat)

    verbose = length(tels) > 10_000
    verbose && print("dots out of 20: ")
    todo, done, pctg = length(tels), 0, 0
    locmat = zeros(BEAST.scalartype(biop, trefs, brefs), BEAST.numfunctions(trefs), numfunctions(brefs))
    for (p,cell) in enumerate(tels)
        P = ta2g[p]
        q = bg2a[P]
        q == 0 && continue

        qr = BEAST.quadrule(biop, trefs, brefs, cell, qd, quadstrat)
        fill!(locmat, 0)
        BEAST.cellinteractions_matched!(locmat, biop, trefs, brefs, cell, qr,p)

        for i in 1 : size(locmat, 1), j in 1 : size(locmat, 2)
            for (m,a) in tad[p,i], (n,b) in bad[q,j]
                store(a * locmat[i,j] * b, m, n)
        
        end end

        new_pctg = round(Int, (done += 1) / todo * 100)
        verbose && new_pctg > pctg + 4 && (print("."); pctg = new_pctg)
    end
end

function cellinteractions(biop::ConductivityTDOp, trefs::U, brefs::V, cell, qr, p) where {U<:RefSpace{T},V<:RefSpace{T}} where {T}

    num_tshs = length(qr[1][3])
    num_bshs = length(qr[1][4])

    zlocal = zeros(T, num_tshs, num_bshs)
    for (i,q) in enumerate(qr)

        w, mp, tvals, bvals = q[1], q[2], q[3], q[4]
        j = w * BEAST.jacobian(mp)
        kernel = BEAST.kernelvals(biop, mp, p, i)

        for m in 1 : num_tshs
            tval = tvals[m]

            for n in 1 : num_bshs
                bval = bvals[n]

                igd = BEAST.integrand(biop, kernel, mp, tval, bval)
                zlocal[m,n] += j * igd

            end
        end
    end

    return zlocal
end

function cellinteractions_matched!(zlocal, biop::ConductivityTDOp, trefs, brefs, cell, qr, p)

    num_tshs = length(qr[1][3])
    num_bshs = length(qr[1][4])

    # zlocal = zeros(Float64, num_tshs, num_bshs)
    for (i,q) in enumerate(qr)

        w, mp, tvals, bvals = q[1], q[2], q[3], q[4]
        j = w * BEAST.jacobian(mp)
        kernel = BEAST.kernelvals(biop, mp, p,i)
        
        for n in 1 : num_bshs
            bval = bvals[n]
            for m in 1 : num_tshs
                tval = tvals[m]

                igd = BEAST.integrand(biop, kernel, mp, tval, bval)
                zlocal[m,n] += j * igd
            end
        end
    end

    return zlocal
end

function assemble!(field::ConductivityTDFunc, tfs::BEAST.Space, store;
    quadstrat=BEAST.defaultquadstrat(field, tfs))

    tels, tad = BEAST.assemblydata(tfs)

    trefs = BEAST.refspace(tfs)
    qd = BEAST.quaddata(field, trefs, tels, quadstrat)

    for (t, tcell) in enumerate(tels)

        # compute the testing with the reference elements
        qr = BEAST.quadrule(field, trefs, t, tcell, qd, quadstrat)
        blocal = BEAST.celltestvalues(trefs, t, tcell, field, qr)

        for i in 1 : BEAST.numfunctions(trefs)
            for (m,a) in tad[t,i]
                store(a*blocal[i], m)
            end
        end

    end

end

function celltestvalues(tshs::BEAST.RefSpace{T, NF}, t, tcell, field::ConductivityTDFunc, qr) where {T, NF}

    num_tshs = numfunctions(tshs)
    interactions = zeros(Complex{T}, num_tshs)

    num_oqp = length(qr)

    for p in 1 : num_oqp
        mp = qr[p].point

        dx =qr[p].weight

        fval = field(t,p,mp)
        tvals = qr[p].value

        for m in 1 : num_tshs
            tval = tvals[m]

            igd = BEAST.integrand(field, tval, fval)
            interactions[m] += igd * dx
        end
    end

    return interactions
end