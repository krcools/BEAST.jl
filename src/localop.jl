

using CollisionDetection

mutable struct SingleQuadStrategy{T}
    coords::Vector{T}
    weights::Vector{T}
end

abstract type LocalOperator <: Operator end


function assemble!(biop::LocalOperator, tfs::Space, bfs::Space, store)

    if geometry(tfs) == geometry(bfs)
        return assemble_local_matched!(biop, tfs, bfs, store)
    end

    if CompScienceMeshes.refines(geometry(tfs), geometry(bfs))
        return assemble_local_refines!(biop, tfs, bfs, store)
    end

    # tels, tad = assemblydata(tfs)
    # bels, bad = assemblydata(bfs)
    #
    # length(tels) != numcells(geometry(tfs)) && return assemble_local_mixed!(biop, tfs, bfs, store)
    # length(bels) != numcells(geometry(bfs)) && return assemble_local_mixed!(biop, tfs, bfs, store)

    return assemble_local_mixed!(biop, tfs, bfs, store)
end

function assemble_local_matched!(biop::LocalOperator, tfs::Space, bfs::Space, store)

    tels, tad, ta2g = assemblydata(tfs)
    bels, bad, ba2g = assemblydata(bfs)

    bg2a = zeros(Int, length(geometry(bfs)))
    for (i,j) in enumerate(ba2g) bg2a[j] = i end

    trefs = refspace(tfs)
    brefs = refspace(bfs)

    qd = quaddata(biop, trefs, brefs, tels, bels)
    for (p,cell) in enumerate(tels)
        P = ta2g[p]
        q = bg2a[P]
        q == 0 && continue

        qr = quadrule(biop, trefs, brefs, cell, qd)
        locmat = cellinteractions(biop, trefs, brefs, cell, qr)

        for i in 1 : size(locmat, 1), j in 1 : size(locmat, 2)
            for (m,a) in tad[p,i], (n,b) in bad[q,j]
                store(a * locmat[i,j] * b, m, n)
end end end end


function assemble_local_refines!(biop::LocalOperator, tfs::Space, bfs::Space, store)

    println("Using 'refines' algorithm for local assembly:")

    # tol = sqrt(eps(Float64))
    tgeo = geometry(tfs)
    bgeo = geometry(bfs)
    @assert CompScienceMeshes.refines(tgeo, bgeo)

    trefs = refspace(tfs)
    brefs = refspace(bfs)

    tels, tad, ta2g = assemblydata(tfs)
    bels, bad, ba2g = assemblydata(bfs)

    bg2a = zeros(Int, length(geometry(bfs)))
    for (i,j) in enumerate(ba2g) bg2a[j] = i end

    qd = quaddata(biop, trefs, brefs, tels, bels)

    # store the bcells in an octree
    # tree = elementstree(bels)

    print("dots out of 10: ")
    todo, done, pctg = length(tels), 0, 0
    for (p,tcell) in enumerate(tels)

        P = ta2g[p]
        Q = CompScienceMeshes.parent(tgeo, P)
        q = bg2a[Q]

        # tc, ts = boundingbox(tcell.vertices)
        # pred = (c,s) -> boxesoverlap(c,s,tc,ts)

        # for box in boxes(tree, pred)
        #     for q in box
        bcell = bels[q]
        @assert overlap(tcell, bcell)

        # if overlap(tcell, bcell)

        isct = intersection(tcell, bcell)
        for cell in isct
            # volume(cell) < tol && continue

            P = restrict(brefs, bcell, cell)
            Q = restrict(trefs, tcell, cell)

            qr = quadrule(biop, trefs, brefs, cell, qd)
            zlocal = cellinteractions(biop, trefs, brefs, cell, qr)
            zlocal = Q * zlocal * P'

            for i in 1 : numfunctions(trefs)
                for j in 1 : numfunctions(brefs)
                    for (m,a) in tad[p,i]
                        for (n,b) in bad[q,j]
                            store(a * zlocal[i,j] * b, m, n)
                        end # next basis function this cell supports
                    end # next test function this cell supports
                end # next refshape on basis side
            end # next refshape on test side

        end # next cell in intersection
        # end # if overlap
        #     end # next cell in the basis geometry
        # end # next box in the octree

        done += 1
        new_pctg = round(Int, done / todo * 100)
        if new_pctg > pctg + 9
            print(".")
            pctg = new_pctg
        end
    end # next cell in the test geometry

    println("")

end

function assemble_local_matched!(biop::LocalOperator, tfs::subdBasis, bfs::subdBasis, store)

    tels, tad = assemblydata(tfs)
    bels, bad = assemblydata(bfs)

    trefs = refspace(tfs)
    brefs = refspace(bfs)

    qd = quaddata(biop, trefs, brefs, tels, bels)
    for (p,cell) in enumerate(tels)

        qr = quadrule(biop, trefs, brefs, cell, qd)
        locmat = cellinteractions(biop, trefs, brefs, cell, qr)

        for i in 1 : size(locmat, 1), j in 1 : size(locmat, 2)
            for (m,a) in tad[p][i], (n,b) in bad[p][j]
                store(a * locmat[i,j] * b, m, n)
end end end end


function elementstree(elements)

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

    return Octree(points, radii)
end


"""
    assemble_local_mixed(biop::LocalOperator, tfs, bfs)

For use when basis and test functions are defined on different meshes
"""
function assemble_local_mixed!(biop::LocalOperator, tfs::Space, bfs::Space, store)

    tol = sqrt(eps(Float64))

    trefs = refspace(tfs)
    brefs = refspace(bfs)

    tels, tad = assemblydata(tfs)
    bels, bad = assemblydata(bfs)

    qd = quaddata(biop, trefs, brefs, tels, bels)

    # store the bcells in an octree
    tree = elementstree(bels)

    print("dots out of 10: ")
    todo, done, pctg = length(tels), 0, 0
    for (p,tcell) in enumerate(tels)

        tc, ts = boundingbox(tcell.vertices)
        pred = (c,s) -> boxesoverlap(c,s,tc,ts)

        for box in boxes(tree, pred)
            for q in box
                bcell = bels[q]

                if overlap(tcell, bcell)

                    isct = intersection(tcell, bcell)
                    for cell in isct
                        volume(cell) < tol && continue

                        P = restrict(brefs, bcell, cell)
                        Q = restrict(trefs, tcell, cell)

                        qr = quadrule(biop, trefs, brefs, cell, qd)
                        zlocal = cellinteractions(biop, trefs, brefs, cell, qr)
                        zlocal = Q * zlocal * P'

                        for i in 1 : numfunctions(trefs)
                            for j in 1 : numfunctions(brefs)
                                for (m,a) in tad[p,i]
                                    for (n,b) in bad[q,j]
                                        store(a * zlocal[i,j] * b, m, n)
                                    end # next basis function this cell supports
                                end # next test function this cell supports
                            end # next refshape on basis side
                        end # next refshape on test side

                    end # next cell in intersection
                end # if overlap
            end # next cell in the basis geometry
        end # next box in the octree

        done += 1
        new_pctg = round(Int, done / todo * 100)
        if new_pctg > pctg + 9
            print(".")
            pctg = new_pctg
        end
    end # next cell in the test geometry

    println("")
end

function cellinteractions(biop, trefs, brefs, cell, qr)

    # num_tshs = numfunctions(trefs)
    num_tshs = length(qr[1][3])
    # num_bshs = numfunctions(brefs)
    num_bshs = length(qr[1][4])

    zlocal = zeros(Float64, num_tshs, num_bshs)
    for q in qr

        w, mp, tvals, bvals = q[1], q[2], q[3], q[4]
        j = w * jacobian(mp)

        kernel = kernelvals(biop, mp)

        # for m in 1 : num_tshs
        for m in 1 : length(tvals)

            tval = tvals[m]
            # for n in 1 : num_bshs
            for n in 1 : length(bvals)
                bval = bvals[n]

                igd = integrand(biop, kernel, mp, tval, bval)
                zlocal[m,n] += j * igd

            end
        end
    end

    return zlocal
end

function allocatestorage(operator::LocalOperator, test_functions, trial_functions,
    ::Type{Val{:bandedstorage}},
    ::Type{BEAST.LongDelays{:ignore}})

    T = promote_type(
        scalartype(operator)       ,
        scalartype(test_functions) ,
        scalartype(trial_functions),
    )
    Z = sparse(SharedArray{T}(
        numfunctions(test_functions)  ,
        numfunctions(trial_functions),
    ))
    store(v,m,n) = (Z[m,n] += v)
    return Z, store
end
