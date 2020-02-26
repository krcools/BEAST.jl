using LinearAlgebra


function addf!(fn::Vector{<:Shape}, x::Vector, space::Space, idcs::Vector{Int})
    for (m,bf) in enumerate(space.fns)
        for sh in bf
            cellid = idcs[sh.cellid]
            BEAST.add!(fn, cellid, sh.refid, sh.coeff * x[m])
        end
    end
end

function Base.in(mesh::CompScienceMeshes.AbstractMesh)
    cells_mesh = sort.(mesh)
    function f(cell)
        sort(cell) in cells_mesh
    end
end



function builddual2form(support, port, dirichlet, prt_fluxes)

    faces = CompScienceMeshes.skeleton_fast(support,2)
    edges = CompScienceMeshes.skeleton_fast(support,1)
    verts = CompScienceMeshes.skeleton_fast(support,0)
    @assert length(verts) - length(edges) + length(faces) - length(support) == 1
    bndry = boundary(support)

    cells_bndry = sort.(bndry)
    dirbnd = submesh(dirichlet) do face
        sort(face) in cells_bndry ? true : false
    end

    bnd_dirbnd = boundary(dirbnd)
    edges_dirbnd = CompScienceMeshes.skeleton_fast(dirbnd,1)
    cells_bnd_dirbnd = sort.(bnd_dirbnd)
    int_edges_dirbnd = submesh(edges_dirbnd) do edge
        sort(edge) in cells(bnd_dirbnd) ? false : true
    end

    num_faces_on_port = 0
    num_faces_on_dirc = 0

    cells_port = sort.(port)
    cells_dirbnd = sort.(dirbnd)
    int_faces = submesh(faces) do face
        sort(face) in cells_port && return false
        sort(face) in cells_dirbnd && return true
        sort(face) in cells_bndry && return false
        return true
    end

    bnd_edges = CompScienceMeshes.skeleton_fast(bndry,1)
    prt_edges = CompScienceMeshes.skeleton_fast(port,1)

    cells_int_edges_dirbnd = sort.(int_edges_dirbnd)
    cells_bnd_edges = sort.(bnd_edges)
    cells_prt_edges = sort.(prt_edges)

    int_edges = submesh(edges) do edge
        sort(edge) in cells_int_edges_dirbnd && return true
        sort(edge) in cells_bnd_edges && return false
        sort(edge) in cells_prt_edges && return false
        return true
    end

    RT_int = nedelecd3d(support, int_faces)
    RT_prt = nedelecd3d(support, port)
    L0_int = nedelecc3d(support, int_edges)

    Id = BEAST.Identity()
    div_RT_int = divergence(RT_int)
    div_RT_prt = divergence(RT_prt)
    D = assemble(Id, div_RT_int, div_RT_int)
    Q = assemble(Id, div_RT_int, div_RT_prt)
    d = -Q * prt_fluxes

    curl_L0_int = curl(L0_int)
    div_curl_L0_int = divergence(curl_L0_int)
    ZZ = real(assemble(Id, div_curl_L0_int, div_curl_L0_int))
    @assert isapprox(norm(ZZ), 0.0, atol=1e-8)
    C = assemble(Id, curl_L0_int, RT_int)
    c = real(assemble(Id, curl_L0_int, RT_prt)) * prt_fluxes

    T = eltype(D)
    nz = length(c)
    QQ = [D C'; C zeros(T,nz,nz)]
    qq = [d;c]
    x = (QQ \ qq)[1:end-nz]

    if !isapprox(C*x, c, atol=1e-8) || !isapprox(D*x, d, atol=1e-6)
        @show norm(D*x-d)
        @show norm(C*x-c)
        @show rank(C)
        @show size(C,1)
        error("error")
    end

    return RT_int, RT_prt, x, prt_fluxes

end


function dual2forms_init(Tetrs)

    tetrs = barycentric_refinement(Tetrs)
    v2t, v2n = CompScienceMeshes.vertextocellmap(tetrs)
    bnd = boundary(tetrs)

    return tetrs, bnd, v2t, v2n
end


function dual2forms(Tetrs, Edges, Dir)
    tetrs, bnd, v2t, v2n = dual2forms_init(Tetrs)
    dual2forms_body(Tetrs, Edges, Dir, tetrs, bnd, v2t, v2n)
end

function dual2forms_body(Tetrs, Edges, Dir, tetrs, bnd, v2t, v2n)

    T = coordtype(Tetrs)
    bfs = Vector{Vector{Shape{T}}}(undef, numcells(Edges))
    pos = Vector{vertextype(Edges)}(undef, numcells(Edges))

    gpred = CompScienceMeshes.overlap_gpredicate(Dir)
    dirichlet = submesh(bnd) do face
        gpred(chart(bnd,face))
    end

    for (F,Edge) in enumerate(cells(Edges))

        println("Constructing dual 2-forms: $F out of $(length(Edges)).")
        bfs[F] = Vector{Shape{T}}()

        pos[F] = cartesian(center(chart(Edges,Edge)))
        port_vertex_idx = argmin(norm.(vertices(tetrs) .- Ref(pos[F])))

        ptch_idcs1 = v2t[Edge[1],1:v2n[Edge[1]]]
        ptch_idcs2 = v2t[Edge[2],1:v2n[Edge[2]]]

        patch1 = Mesh(vertices(tetrs), cells(tetrs)[ptch_idcs1])
        patch2 = Mesh(vertices(tetrs), cells(tetrs)[ptch_idcs2])
        patch = CompScienceMeshes.union(patch1, patch2)

        patch_bnd = boundary(patch1)

        bnd_patch1 = boundary(patch1)
        bnd_patch2 = boundary(patch2)
        port = submesh(in(bnd_patch2), bnd_patch1)

        total_area = sum(volume(chart(port, fc)) for fc in port)
        # prt_fluxes = ones(T, numcells(port)) / numcells(port)
        prt_fluxes = zeros(length(port))
        tgt = vertices(Edges)[Edge[1]] - vertices(Edges)[Edge[2]]
        for (i,face) in enumerate(cells(port))
            chrt = chart(port, face)
            sgn = sign(dot(normal(chrt), tgt))
            # prt_fluxes[i] *= sgn
            prt_fluxes[i] = sgn * volume(chrt) / total_area
        end

        # RT_int, RT_prt, x_int, x_prt = builddual2form(patch, port, dirichlet, prt_fluxes)
        # ptch_idcs = [ptch_idcs1; ptch_idcs2]
        # addf!(bfs[F], x_int, RT_int, ptch_idcs)
        # addf!(bfs[F], x_prt, RT_prt, ptch_idcs)

        RT1_int, RT1_prt, x1_int, x_prt = builddual2form(
            patch1, port, dirichlet, +prt_fluxes)
        RT2_int, RT2_prt, x2_int, x_prt = builddual2form(
            patch2, port, dirichlet, prt_fluxes)

        addf!(bfs[F], x1_int, RT1_int, ptch_idcs1)
        addf!(bfs[F], +x_prt, RT1_prt, ptch_idcs1)

        addf!(bfs[F], x2_int, RT2_int, ptch_idcs2)
        addf!(bfs[F], x_prt, RT2_prt, ptch_idcs2)

    end

    NDLCDBasis(tetrs, bfs, pos)
end


function extend_edge_to_face(supp, dirichlet, x_prt, port_edges)

    Id = BEAST.Identity()

    bnd_supp = boundary(supp)
    supp_edges = CompScienceMeshes.skeleton_fast(supp, 1)
    supp_nodes = CompScienceMeshes.skeleton_fast(supp, 0)

    dir_compl = submesh(!in(dirichlet), bnd_supp)
    # dir_compl_edges = submesh(!in(dirichlet), bnd_supp)
    dir_compl_edges = CompScienceMeshes.skeleton_fast(dir_compl, 1)
    dir_compl_nodes = CompScienceMeshes.skeleton_fast(dir_compl, 0)

    int_edges = submesh(!in(dir_compl_edges), supp_edges)
    int_nodes = submesh(!in(dir_compl_nodes), supp_nodes)

    Nd_prt = BEAST.nedelec(supp, port_edges)
    Nd_int = BEAST.nedelec(supp, int_edges)
    Lg_int = BEAST.lagrangec0d1(supp, int_nodes)

    curl_Nd_prt = divergence(n × Nd_prt)
    curl_Nd_int = divergence(n × Nd_int)
    grad_Lg_int = n × curl(Lg_int)

    A = assemble(Id, curl_Nd_int, curl_Nd_int)
    B = assemble(Id, grad_Lg_int, Nd_int)

    a = -assemble(Id, curl_Nd_int, curl_Nd_prt) * x_prt
    b = -assemble(Id, grad_Lg_int, Nd_prt) * x_prt

    Z = zeros(eltype(b), length(b), length(b))
    u = [A B'; B Z] \ [a;b]
    x_int = u[1:end-length(b)]

    return x_int, int_edges, Nd_int
end



function extend_face_to_tetr(supp, dirichlet, x_prt, port_edges)

    Id = BEAST.Identity()

    bnd_supp = boundary(supp)
    supp_edges = CompScienceMeshes.skeleton_fast(supp, 1)
    supp_nodes = CompScienceMeshes.skeleton_fast(supp, 0)

    dir_compl = submesh(!in(dirichlet), bnd_supp)
    dir_compl_edges = CompScienceMeshes.skeleton_fast(dir_compl, 1)
    dir_compl_nodes = CompScienceMeshes.skeleton_fast(dir_compl, 0)

    int_edges = submesh(!in(dir_compl_edges), supp_edges)
    int_nodes = submesh(!in(dir_compl_nodes), supp_nodes)

    Nd_prt = BEAST.nedelecc3d(supp, port_edges)
    Nd_int = BEAST.nedelecc3d(supp, int_edges)
    Lg_int = BEAST.lagrangec0d1(supp, int_nodes)

    curl_Nd_prt = curl(Nd_prt)
    curl_Nd_int = curl(Nd_int)
    grad_Lg_int = BEAST.gradient(Lg_int)

    A = assemble(Id, curl_Nd_int, curl_Nd_int)
    B = assemble(Id, grad_Lg_int, Nd_int)

    a = -assemble(Id, curl_Nd_int, curl_Nd_prt) * x_prt
    b = -assemble(Id, grad_Lg_int, Nd_prt) * x_prt

    Z = zeros(eltype(b), length(b), length(b))
    u = [A B'; B Z] \ [a;b]
    x_int = u[1:end-length(b)]

    return x_int, int_edges, Nd_int
end

function dual1forms(Tetrs, Faces, Dir)
    tetrs, bnd, dir, v2t, v2n = dual1forms_init(Tetrs, Dir)
    dual1forms_body(Faces, tetrs, bnd, dir, v2t, v2n)
end

function dual1forms_init(Tetrs, Dir)
    tetrs = barycentric_refinement(Tetrs)
    v2t, v2n = CompScienceMeshes.vertextocellmap(tetrs)
    bnd = boundary(tetrs)
    gpred = CompScienceMeshes.overlap_gpredicate(Dir)
    dir = submesh(face -> gpred(chart(bnd,face)), bnd)
    return tetrs, bnd, dir, v2t, v2n
end

function dual1forms_body(Faces, tetrs, bnd, dir, v2t, v2n)

    T = coordtype(tetrs)
    bfs = Vector{Vector{Shape{T}}}(undef, length(Faces))
    pos = Vector{vertextype(Faces)}(undef, length(Faces))

    for (F,Face) in enumerate(Faces)

        println("Constructing dual 1-forms: $F out of $(length(Faces)).")

        idcs1 = v2t[Face[1],1:v2n[Face[1]]]
        idcs2 = v2t[Face[2],1:v2n[Face[2]]]
        idcs3 = v2t[Face[3],1:v2n[Face[3]]]

        supp1 = tetrs.mesh[idcs1]
        supp2 = tetrs.mesh[idcs2]
        supp3 = tetrs.mesh[idcs3]

        bnd_supp1 = boundary(supp1)
        bnd_supp2 = boundary(supp2)
        bnd_supp3 = boundary(supp3)

        dir1_faces = submesh(in(dir), bnd_supp1)
        dir2_faces = submesh(in(dir), bnd_supp2)
        dir3_faces = submesh(in(dir), bnd_supp3)

        bnd_dir1 = boundary(dir1_faces)
        bnd_dir2 = boundary(dir2_faces)
        bnd_dir3 = boundary(dir3_faces)

        supp23 = submesh(in(bnd_supp2), bnd_supp3)
        supp31 = submesh(in(bnd_supp3), bnd_supp1)
        supp12 = submesh(in(bnd_supp1), bnd_supp2)

        dir23_edges = submesh(in(bnd_dir2), bnd_dir3)
        dir31_edges = submesh(in(bnd_dir3), bnd_dir1)
        dir12_edges = submesh(in(bnd_dir1), bnd_dir2)

        port_edges = boundary(supp23)
        port_edges = submesh(in(boundary(supp31)), port_edges)
        port_edges = submesh(in(boundary(supp12)), port_edges)
        @assert 1 ≤ length(port_edges) ≤ 2

        # Step 1: set port flux and extend to dual faces
        x0 = zeros(length(port_edges))
        total_lgt = sum(volume(chart(port_edges, edge)) for edge in port_edges)
        for (i,edge) in enumerate(port_edges)
            tgt = tangents(center(chart(port_edges, edge)),1)
            lgt = volume(chart(port_edges, edge))
            sgn = sign(dot(normal(chart(Faces, Face)), tgt))
            x0[i] = sgn * lgt / total_lgt
        end

        x23, supp23_int_edges = extend_edge_to_face(supp23, dir23_edges, x0, port_edges)
        x31, supp31_int_edges = extend_edge_to_face(supp31, dir31_edges, x0, port_edges)
        x12, supp12_int_edges = extend_edge_to_face(supp12, dir12_edges, x0, port_edges)

        port1_edges = CompScienceMeshes.union(port_edges, supp31_int_edges)
        port1_edges = CompScienceMeshes.union(port1_edges, supp12_int_edges)
        port2_edges = CompScienceMeshes.union(port_edges, supp12_int_edges)
        port2_edges = CompScienceMeshes.union(port2_edges, supp23_int_edges)
        port3_edges = CompScienceMeshes.union(port_edges, supp23_int_edges)
        port3_edges = CompScienceMeshes.union(port3_edges, supp31_int_edges)

        x1_prt = [x0; x31; x12]
        x2_prt = [x0; x12; x23]
        x3_prt = [x0; x23; x31]

        Nd1_prt = BEAST.nedelecc3d(supp1, port1_edges)
        Nd2_prt = BEAST.nedelecc3d(supp2, port2_edges)
        Nd3_prt = BEAST.nedelecc3d(supp3, port3_edges)

        x1_int, _, Nd1_int = extend_face_to_tetr(supp1, dir1_faces, x1_prt, port1_edges)
        x2_int, _, Nd2_int = extend_face_to_tetr(supp2, dir2_faces, x2_prt, port2_edges)
        x3_int, _, Nd3_int = extend_face_to_tetr(supp3, dir3_faces, x3_prt, port3_edges)

        # inject in the global space
        fn = BEAST.Shape{Float64}[]
        addf!(fn, x1_prt, Nd1_prt, idcs1)
        addf!(fn, x1_int, Nd1_int, idcs1)

        addf!(fn, x2_prt, Nd2_prt, idcs2)
        addf!(fn, x2_int, Nd2_int, idcs2)

        addf!(fn, x3_prt, Nd3_prt, idcs3)
        addf!(fn, x3_int, Nd3_int, idcs3)

        pos[F] = cartesian(CompScienceMeshes.center(chart(Faces, Face)))
        bfs[F] = fn
        # space = BEAST.NDLCCBasis(tetrs, [fn], [pos])
    end

    NDLCCBasis(tetrs, bfs, pos)
end
