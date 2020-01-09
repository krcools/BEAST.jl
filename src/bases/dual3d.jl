using LinearAlgebra

function builddual2form(support, port, dirichlet, prt_fluxes)

    println()

    faces = skeleton(support,2)
    edges = skeleton(support,1)
    verts = skeleton(support,0)
    @assert length(verts) - length(edges) + length(faces) - length(support) == 1
    bndry = boundary(support)
    @show numcells(faces)

    cells_bndry = [sort(c) for c in cells(bndry)]
    dirbnd = submesh(dirichlet) do face
        sort(face) in cells(bndry) ? true : false
    end
    @show numcells(support)
    @show numcells(dirbnd)
    # @assert numcells(dirbnd) ≤ 4

    bnd_dirbnd = boundary(dirbnd)
    edges_dirbnd = skeleton(dirbnd,1)
    cells_bnd_dirbnd = [sort(c) for c in cells(bnd_dirbnd)]
    int_edges_dirbnd = submesh(edges_dirbnd) do edge
        sort(edge) in cells(bnd_dirbnd) ? false : true
    end
    @show numcells(int_edges_dirbnd)
    # @assert numcells(int_edges_dirbnd) ≤ 2

    int_pred = interior_tpredicate(support)
    num_faces_on_port = 0
    num_faces_on_dirc = 0

    @assert numcells(submesh(!int_pred, faces)) == numcells(boundary(support))


    # for edge in cells(faces)
    #     println(edge)
    # end
    # println()
    # for edge in cells(dirbnd)
    #     println(edge)
    # end
    # println()
    # for edge in cells(port)
    #     println(edge)
    # end

    cells_port = [sort(c) for c in cells(port)]
    cells_dirbnd = [sort(c) for c in cells(dirbnd)]
    int_faces = submesh(faces) do face
        sort(face) in cells_port && return false
        sort(face) in cells_dirbnd && return true
        int_pred(face) ? true : false
    end
    # println()
    # for edge in cells(int_faces)
    #     println(edge)
    # end
    @show numcells(int_faces)
    # @show num_faces_on_port
    # @show num_faces_on_dirc

    bnd_edges = skeleton(bndry,1)
    prt_edges = skeleton(port,1)

    cells_int_edges_dirbnd = [sort(c) for c in cells(int_edges_dirbnd)]
    cells_bnd_edges = [sort(c) for c in cells(bnd_edges)]
    cells_prt_edges = [sort(c) for c in cells(prt_edges)]

    @show length(cells_int_edges_dirbnd)
    @show length(cells_bnd_edges)
    @show length(cells_prt_edges)
    int_edges = submesh(edges) do edge
        sort(edge) in cells_int_edges_dirbnd && return true
        sort(edge) in cells_bnd_edges && return false
        sort(edge) in cells_prt_edges && return false
        return true
    end
    @show numcells(int_edges)

    RT_int = nedelecd3d(support, int_faces)
    RT_prt = nedelecd3d(support, port)
    # vertex_list = [c[1] for c in cells(int_edges)]
    # L0_int = lagrangec0d1(suppport, vertex_list, Val{4})
    L0_int = nedelecc3d(support, int_edges)

    Id = BEAST.Identity()
    D = assemble(Id, divergence(RT_int), divergence(RT_int))
    Q = assemble(Id, divergence(RT_int), divergence(RT_prt))
    d = -Q * prt_fluxes

    @show numfunctions(L0_int)
    # @assert numfunctions(L0_int) in [1,2]
    curl_L0_int = curl(L0_int)
    div_curl_L0_int = divergence(curl_L0_int)
    ZZ = real(assemble(Id, div_curl_L0_int, div_curl_L0_int))
    @assert isapprox(norm(ZZ), 0.0, atol=1e-8)
    C = assemble(Id, curl_L0_int, RT_int)
    c = real(assemble(Id, curl_L0_int, RT_prt)) * prt_fluxes

    x1 = pinv(D) * d
    N = nullspace(D)
    @show size(N)
    @show rank(C)
    @assert size(N,2) == rank(C)
    @show size(C)
    # @assert rank(C) == size(C,1)
    p = (C*N) \ (c - C*x1)
    x = x1 + N*p

    # D*x ≈ d atol=1e-8
    if !isapprox(C*x, c, atol=1e-8) || !isapprox(D*x, d, atol=1e-6)
        @show norm(D*x-d)
        @show norm(C*x-c)
        @show rank(C)
        @show size(C,1)
        error("error")
    end

    if rank(C) != size(C,1)
        @show rank(C)
        @show size(C)
    end

    return RT_int, RT_prt, x, prt_fluxes

end


function dual2forms(Tetrs, Edges)

    tetrs = barycentric_refinement(Tetrs)
    # Bndry = boundary(Tetrs)

    T = coordtype(Tetrs)
    bfs = Vector{Vector{Shape{T}}}(undef, numcells(Edges))
    pos = Vector{vertextype(Edges)}(undef, numcells(Edges))
    dirichlet = boundary(tetrs)
    for (F,Edge) in enumerate(cells(Edges))
        @show F
        bfs[F] = Vector{Shape{T}}()

        pos[F] = cartesian(center(chart(Edges,Edge)))
        port_vertex_idx = argmin(norm.(vertices(tetrs) .- Ref(pos[F])))

        # Build the dual support
        ptch_idcs1 = [i for (i,tetr) in enumerate(cells(tetrs)) if Edge[1] in tetr]
        ptch_idcs2 = [i for (i,tetr) in enumerate(cells(tetrs)) if Edge[2] in tetr]
        patch1 = Mesh(vertices(tetrs), cells(tetrs)[ptch_idcs1])
        patch2 = Mesh(vertices(tetrs), cells(tetrs)[ptch_idcs2])
        patch_bnd = boundary(patch1)
        # @assert CompScienceMeshes.isoriented(patch_bnd)
        port = Mesh(vertices(tetrs), filter(c->port_vertex_idx in c, cells(patch_bnd)))
        # @assert CompScienceMeshes.isoriented(port)
        patch = CompScienceMeshes.union(patch1, patch2)
        @show numcells(patch_bnd)
        @show numcells(patch)
        @show numcells(port)
        # @assert numcells(skeleton(patch,1))+2 == numcells(skeleton(patch1,1)) + numcells(skeleton(patch2,1))

        # @assert numcells(patch) >= 6
        # @assert numcells(port) == 2
        prt_fluxes = ones(T, numcells(port)) / numcells(port)
        tgt = vertices(Edges)[Edge[1]] - vertices(Edges)[Edge[2]]
        for (i,face) in enumerate(cells(port))
            chrt = chart(port, face)
            prt_fluxes[i] *= sign(dot(normal(chrt), tgt))
        end
        RT_int, RT_prt, x_int, x_prt = builddual2form(patch, port, dirichlet, prt_fluxes)

        ptch_idcs = vcat(ptch_idcs1, ptch_idcs2)
        for (m,bf) in enumerate(RT_int.fns)
            for sh in bf
                cellid = ptch_idcs[sh.cellid]
                BEAST.add!(bfs[F], cellid, sh.refid, sh.coeff * x_int[m])
            end
        end

        for (m,bf) in enumerate(RT_prt.fns)
            for sh in bf
                cellid = ptch_idcs[sh.cellid]
                BEAST.add!(bfs[F],cellid, sh.refid, sh.coeff * x_prt[m])
            end
        end

        # F == 25 && error("stop")
    end

    NDLCDBasis(tetrs, bfs, pos)
end
