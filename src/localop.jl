

using CollisionDetection


abstract type LocalOperator <: Operator end


function allocatestorage(op::LocalOperator, test_functions, trial_functions,
    storage_trait::Type{Val{:bandedstorage}})

    T = scalartype(op, test_functions, trial_functions)

    M = Int[]
    N = Int[]
    V = T[]

    function store(v,m,n)
        push!(M,m)
        push!(N,n)
        push!(V,v)
    end

    function freeze()
        nrows = numfunctions(test_functions)
        ncols = numfunctions(trial_functions)
        return sparse(M,N,V, nrows, ncols)
    end

    return freeze, store
end

function allocatestorage(op::LocalOperator, testfunctions, trialfunctions,
    storage_trait::Type{Val{:sparsedicts}})

    T = scalartype(op, testfunctions, trialfunctions)

    m = numfunctions(testfunctions)
    n = numfunctions(trialfunctions)
    Z = ExtendableSparseMatrix(T,m,n)

    store(v,m,n) = (Z[m,n] += v)
    freeze() = SparseArrays.SparseMatrixCSC(Z)

    return freeze, store
end

function allocatestorage(op::LocalOperator, test_functions, trial_functions,
    storage_trait::Type{Val{:densestorage}})

    T = scalartype(op, test_functions, trial_functions)

    Z = zeros(T, numfunctions(test_functions), numfunctions(trial_functions))
    store(v,m,n) = (Z[m,n] += v)
    freeze() = Z

    return freeze, store
end

function assemble!(op::LocalOperator, tfs::Space, bfs::Space, store,
        threading::Type{Threading{:cellcoloring}};
        quadstrat=defaultquadstrat,
        kwargs...)

        assemble!(op, tfs, bfs, store, Threading{:multi}; quadstrat, kwargs...)
end

function assemble!(biop::LocalOperator, tfs::Space, bfs::Space, store,
        threading::Type{Threading{:multi}};
        quadstrat=defaultquadstrat,
        kwargs...)

        quadstrat = quadstrat(biop, tfs, bfs)

    numfunctions(tfs) == 0 && return
    numfunctions(bfs) == 0 && return
 
    if geometry(tfs) == geometry(bfs)
        return assemble_local_matched!(biop, tfs, bfs, store; quadstrat)
    end

    if CompScienceMeshes.refines(geometry(tfs), geometry(bfs))
        return assemble_local_refines!(biop, tfs, bfs, store; quadstrat)
    end

    return assemble_local_mixed!(biop, tfs, bfs, store; quadstrat)
end

function _scatterlocal!(store, locmat, tad, bad, p, q)
    for i in axes(locmat, 1), j in axes(locmat, 2)
        for (m, a) in tad[p, i], (n, b) in bad[q, j]
            store(a * locmat[i, j] * b, m, n)
        end
    end
    return nothing
end

function _scatterlocal!(store, locmat, tad::AbstractVector, bad::AbstractVector, p, q)
    for i in axes(locmat, 1), j in axes(locmat, 2)
        for (m, a) in tad[p][i], (n, b) in bad[q][j]
            store(a * locmat[i, j] * b, m, n)
        end
    end
    return nothing
end

function _assemble_restricted_interaction!(
    biop,
    trefs,
    brefs,
    tcell,
    bcell,
    cell,
    tad,
    bad,
    p,
    q,
    qd,
    quadstrat,
    store,
    tol,
)
    tol === nothing || volume(cell) < tol && return nothing

    testrestriction = restrict(trefs, tcell, cell)
    trialrestriction = restrict(brefs, bcell, cell)
    qr = quadrule(biop, trefs, brefs, cell, qd, quadstrat)
    zlocal = cellinteractions(biop, trefs, brefs, cell, qr)
    zlocal = testrestriction * zlocal * trialrestriction'
    _scatterlocal!(store, zlocal, tad, bad, p, q)
    return nothing
end

function assemble!(biop::LocalOperator, tfs::Space, bfs::Space, store,
    threading::Type{Threading{:single}};
    quadstrat=defaultquadstrat,
    kwargs...)

    quadstrat = quadstrat(biop, tfs, bfs)

    if geometry(tfs) == geometry(bfs)
        return assemble_local_matched!(biop, tfs, bfs, store; quadstrat)
    end

    if CompScienceMeshes.refines(geometry(tfs), geometry(bfs))
        return assemble_local_refines!(biop, tfs, bfs, store; quadstrat)
    end

    return assemble_local_mixed!(biop, tfs, bfs, store; quadstrat)
end

function assemble_local_matched!(biop::LocalOperator, tfs::Space, bfs::Space, store;
    quadstrat=defaultquadstrat(biop, tfs, bfs))

    tels, tad, ta2g = assemblydata(tfs)
    bels, bad, ba2g = assemblydata(bfs)

    bg2a = zeros(Int, length(geometry(bfs)))
    for (i,j) in enumerate(ba2g) bg2a[j] = i end

    trefs = refspace(tfs)
    brefs = refspace(bfs)

    tgeo = geometry(tfs)
    bgeo = geometry(bfs)

    tdom = domain(chart(tgeo, first(tgeo)))
    bdom = domain(chart(bgeo, first(bgeo)))

    num_trefs = numfunctions(trefs, tdom)
    num_brefs = numfunctions(brefs, bdom)

    qd = quaddata(biop, trefs, brefs, tels, bels, quadstrat)

    verbose = length(tels) > 10_000
    verbose && print("dots out of 20: ")
    todo, done, pctg = length(tels), 0, 0
    locmat = zeros(scalartype(biop, trefs, brefs), num_trefs, num_brefs)
    for (p,cell) in enumerate(tels)
        P = ta2g[p]
        q = bg2a[P]
        q == 0 && continue

        qr = quadrule(biop, trefs, brefs, cell, qd, quadstrat)
        fill!(locmat, 0)
        cellinteractions_matched!(locmat, biop, trefs, brefs, cell, qr)

        _scatterlocal!(store, locmat, tad, bad, p, q)

        new_pctg = round(Int, (done += 1) / todo * 100)
        verbose && new_pctg > pctg + 4 && (print("."); pctg = new_pctg)
    end
end


function assemble_local_refines!(biop::LocalOperator, tfs::Space, bfs::Space, store;
    quadstrat=defaultquadstrat(biop, tfs, bfs))

    # println("Using 'refines' algorithm for local assembly:")

    tgeo = geometry(tfs)
    bgeo = geometry(bfs)
    @assert CompScienceMeshes.refines(tgeo, bgeo)

    trefs = refspace(tfs)
    brefs = refspace(bfs)

    tgeo = geometry(tfs)
    bgeo = geometry(bfs)

    tdom = domain(chart(tgeo, first(tgeo)))
    bdom = domain(chart(bgeo, first(bgeo)))

    tels, tad, ta2g = assemblydata(tfs)
    bels, bad, ba2g = assemblydata(bfs)

    bg2a = zeros(Int, length(geometry(bfs)))
    for (i,j) in enumerate(ba2g) bg2a[j] = i end

    qd = quaddata(biop, trefs, brefs, tels, bels, quadstrat)

    print("dots out of 10: ")
    todo, done, pctg = length(tels), 0, 0
    for (p,tcell) in enumerate(tels)

        P = ta2g[p]
        Q = CompScienceMeshes.parent(tgeo, P)
        q = bg2a[Q]
        q == 0 && continue

        bcell = bels[q]
        @assert overlap(tcell, bcell)

        isct = intersection(tcell, bcell)
        for cell in isct
            _assemble_restricted_interaction!(
                biop,
                trefs,
                brefs,
                tcell,
                bcell,
                cell,
                tad,
                bad,
                p,
                q,
                qd,
                quadstrat,
                store,
                nothing,
            )

        end # next cell in intersection

        done += 1
        new_pctg = round(Int, done / todo * 100)
        if new_pctg > pctg + 9
            print(".")
            pctg = new_pctg
        end
    end # next cell in the test geometry

    println("")

end

function assemble_local_matched!(biop::LocalOperator, tfs::subdBasis, bfs::subdBasis, store;
    quadstrat=defaultquadstrat(biop, tfs, bfs))

    tels, tad = assemblydata(tfs)
    bels, bad = assemblydata(bfs)

    trefs = refspace(tfs)
    brefs = refspace(bfs)

    qd = quaddata(biop, trefs, brefs, tels, bels, quadstrat)
    for (p,cell) in enumerate(tels)

        qr = quadrule(biop, trefs, brefs, cell, qd, quadstrat)
        locmat = cellinteractions(biop, trefs, brefs, cell, qr)

        _scatterlocal!(store, locmat, tad, bad, p, p)
    end
end


function elementstree(elements, expansion_ratio=1.1)

    nverts = dimension(eltype(elements)) + 1
    ncells = length(elements)

    @assert !isempty(elements)
    P = eltype(elements[1].vertices)
    T = coordtype(eltype(elements))

    points = zeros(P, ncells)
    radii = zeros(T, ncells)

    for i in 1 : ncells

        verts = elements[i].vertices

        bary = verts[1]
        for j in 2:length(verts)
            bary += verts[j]
        end

        points[i] = bary / nverts
        for j in 1 : nverts
            radii[i] = max(radii[i], norm(verts[j]-points[i]))
        end
    end

    return Octree(points, radii, T(expansion_ratio))
end

struct MixedAssemblyVisitor{B,C,BR,TR,TA,BA,Q,O,S,T,QS}
    tcell::C
    bels::B
    brefs::BR
    trefs::TR
    tad::TA
    bad::BA
    p::Int
    qd::Q
    biop::O
    store::S
    tol::T
    quadstrat::QS
end

function (visitor::MixedAssemblyVisitor)(q)
    bcell = visitor.bels[q]
    overlap(visitor.tcell, bcell) || return false

    isct = intersection(visitor.tcell, bcell)
    for cell in isct
        _assemble_restricted_interaction!(
            visitor.biop,
            visitor.trefs,
            visitor.brefs,
            visitor.tcell,
            bcell,
            cell,
            visitor.tad,
            visitor.bad,
            visitor.p,
            q,
            visitor.qd,
            visitor.quadstrat,
            visitor.store,
            visitor.tol,
        )
    end
    return false
end

function _assemble_mixed_cell!(
    tree,
    tcell,
    p::Int,
    bels,
    brefs,
    trefs,
    tad,
    bad,
    qd,
    biop,
    store,
    tol,
    quadstrat,
)
    tc, ts = boundingbox(tcell.vertices)
    visitor = MixedAssemblyVisitor(
        tcell,
        bels,
        brefs,
        trefs,
        tad,
        bad,
        p,
        qd,
        biop,
        store,
        tol,
        quadstrat,
    )
    CollisionDetection.foreachsearchtree(visitor, tree, (tc, ts))
    return nothing
end

function _assemble_mixed_cells!(
    tels::AbstractVector{C},
    tree,
    bels,
    brefs,
    trefs,
    tad,
    bad,
    qd,
    biop,
    store,
    tol,
    quadstrat,
) where {C}
    print("dots out of 10: ")
    todo, done, pctg = length(tels), 0, 0
    for (p, tcell) in enumerate(tels)
        _assemble_mixed_cell!(
            tree,
            tcell,
            p,
            bels,
            brefs,
            trefs,
            tad,
            bad,
            qd,
            biop,
            store,
            tol,
            quadstrat,
        )

        done += 1
        new_pctg = round(Int, done / todo * 100)
        if new_pctg > pctg + 9
            print(".")
            pctg = new_pctg
        end
    end
    println("")
    return nothing
end


"""
    assemble_local_mixed(biop::LocalOperator, tfs, bfs)

For use when basis and test functions are defined on different meshes
"""
function assemble_local_mixed!(biop::LocalOperator, tfs::Space{T}, bfs::Space{T}, store;
    quadstrat=defaultquadstrat(biop, tfs, bfs)) where {T}

    tol = sqrt(eps(T))

    trefs = refspace(tfs)
    brefs = refspace(bfs)

    tr = assemblydata(tfs); tr == nothing && return
    br = assemblydata(bfs); br == nothing && return

    tels, tad = tr
    bels, bad = br

    qd = quaddata(biop, trefs, brefs, tels, bels, quadstrat)

    # store the bcells in an octree
    tree = elementstree(bels)

    _assemble_mixed_cells!(
        tels,
        tree,
        bels,
        brefs,
        trefs,
        tad,
        bad,
        qd,
        biop,
        store,
        tol,
        quadstrat,
    )
end


function cellinteractions_matched!(zlocal, biop, trefs, brefs, cell, qr)

    num_tshs = length(qr[1][3])
    num_bshs = length(qr[1][4])

    # zlocal = zeros(Float64, num_tshs, num_bshs)
    for q in qr

        w, mp, tvals, bvals = q[1], q[2], q[3], q[4]
        j = w * jacobian(mp)
        kernel = kernelvals(biop, mp)
        
        for n in 1 : num_bshs
            bval = bvals[n]
            for m in 1 : num_tshs
                tval = tvals[m]

                igd = integrand(biop, kernel, mp, tval, bval)
                zlocal[m,n] += j * igd
            end
        end
    end

    return zlocal
end

function cellinteractions(biop, trefs::U, brefs::V, cell, qr) where {U<:RefSpace{T},V<:RefSpace{T}} where {T}

    num_tshs = length(qr[1][3])
    num_bshs = length(qr[1][4])

    zlocal = zeros(T, num_tshs, num_bshs)
    for q in qr

        w, mp, tvals, bvals = q[1], q[2], q[3], q[4]
        j = w * jacobian(mp)
        kernel = kernelvals(biop, mp)

        for m in 1 : num_tshs
            tval = tvals[m]

            for n in 1 : num_bshs
                bval = bvals[n]

                igd = integrand(biop, kernel, mp, tval, bval)
                zlocal[m,n] += j * igd

            end
        end
    end

    return zlocal
end


@testitem "assemble!: zero sized block" begin
    using CompScienceMeshes

    fn = joinpath(dirname(pathof(BEAST)), "../examples/assets/sphere45.in")
    m1 = readmesh(fn)
    m2 = m1[Int[]]

    X = BEAST.DirectProductSpace([raviartthomas(m) for m in [m1, m2]])
    Id = BEAST.Identity()

    @hilbertspace j[1:2]
    @hilbertspace k[1:2]
    a = Id[k[1],j[1]] + Id[k[2],j[2]]

    A = assemble(a, X, X)
    M = AbstractMatrix(A)
    import BEAST.BlockArrays

    n1 = numfunctions(X[1])
    n2 = numfunctions(X[2])

    @test n2 == 0

    @test BlockArrays.blocksize(M) == (2,2)
    @test BlockArrays.blocksizes(M) == [(n1,n1) (n1,n2); (n2,n1) (n2,n2)]
end
