export localoperator
export localoperator2

using CollisionDetection

type SingleQuadStrategy{T}
    coords::Vector{T}
    weights::Vector{T}
end

abstract type LocalOperator <: Operator end


function assemble!(biop::LocalOperator, tfs::Space, bfs::Space, store)
    if tfs == bfs
        return assemble_local_matched!(biop, tfs, bfs, store)
    else
        return assemble_local_mixed!(biop, tfs, bfs, store)
    end
end

function assemble_local_matched!(biop::LocalOperator, tfs::Space, bfs::Space, store)

    tels, tad = assemblydata(tfs)
    bels, bad = assemblydata(bfs)

    trefs = refspace(tfs)
    brefs = refspace(bfs)

    qd = quaddata(biop, trefs, brefs, tels, bels)

    #for (p,cell) in cellenumeration(geometry(tfs))
    for (p,cell) in enumerate(tels)

        # compute local interactions on reference cell
        qr = quadrule(biop, trefs, brefs, cell, qd)
        locmat = cellinteractions(biop, trefs, brefs, cell, qr)

        # assemble the global matrix
        for i in 1 : size(locmat, 1)
            for j in 1 : size(locmat, 2)

                for (m,a) in tad[p,i]
                    for (n,b) in bad[p,j]
                        store(a * locmat[i,j] * b, m, n)
                    end
                end
            end
        end
    end
end


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

        #points[i] = sum(verts) / nverts
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

    num_tshs = numfunctions(trefs)
    num_bshs = numfunctions(brefs)

    zlocal = zeros(Float64, num_tshs, num_bshs)
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
