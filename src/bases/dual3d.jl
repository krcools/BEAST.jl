using LinearAlgebra


# function addf!(fn::Vector{<:Shape}, x::Vector, space::Space, idcs::Vector{Int})
#     for (m,bf) in enumerate(space.fns)
#         for sh in bf
#             cellid = idcs[sh.cellid]
#             BEAST.add!(fn, cellid, sh.refid, sh.coeff * x[m])
#         end
#     end
# end

function dual2forms_body(Edges, tetrs, bnd, dir, v2t, v2n)

    T = coordtype(tetrs)
    bfs = Vector{Vector{Shape{T}}}(undef, numcells(Edges))
    pos = Vector{vertextype(Edges)}(undef, numcells(Edges))

    Cells = cells(Edges)
    num_threads = Threads.nthreads()
    Threads.@threads for F in eachindex(Cells)
        Edge = Cells[F]

        myid = Threads.threadid()
        myid == 1 && (F % 20 == 0) &&
            println("Constructing dual 2-forms: $(F*num_threads) out of $(length(Edges)).")

        idcs1 = v2t[Edge[1],1:v2n[Edge[1]]]
        idcs2 = v2t[Edge[2],1:v2n[Edge[2]]]

        supp1 = tetrs.mesh[idcs1]
        supp2 = tetrs.mesh[idcs2]

        bnd_supp1 = boundary(supp1)
        bnd_supp2 = boundary(supp2)

        dir1_faces = submesh(in(dir), bnd_supp1)
        dir2_faces = submesh(in(dir), bnd_supp2)

        port_faces = bnd_supp1
        port_faces = submesh(in(bnd_supp2), port_faces)

        x0 = zeros(length(port_faces))
        total_vol = sum(volume(chart(port_faces, fc)) for fc in port_faces)
        tgt = tangents(center(chart(Edges, F)),1)
        for (i,face) in enumerate(port_faces)
            chrt = chart(port_faces, face)
            nrm = normal(chrt)
            vol = volume(chrt)
            sgn = sign(dot(nrm, tgt))
            x0[i] = sgn * vol / total_vol
        end

        RT1_prt = BEAST.raviartthomas(supp1, port_faces)
        RT2_prt = BEAST.raviartthomas(supp2, port_faces)

        x1_int, _, RT1_int = extend_2_form(supp1, dir1_faces, x0, port_faces)
        x2_int, _, RT2_int = extend_2_form(supp2, dir2_faces, x0, port_faces)

        bfs[F] = Vector{Shape{T}}()
        pos[F] = cartesian(center(chart(Edges,F)))
        addf!(bfs[F], x1_int, RT1_int, idcs1)
        addf!(bfs[F], x0, RT1_prt, idcs1)

        addf!(bfs[F], x2_int, RT2_int, idcs2)
        addf!(bfs[F], x0, RT2_prt, idcs2)

    end

    NDLCDBasis(tetrs, bfs, pos)
end


function extend_1_form(supp, dirichlet, x_prt, port_edges)

    Id = BEAST.Identity()

    bnd_supp = boundary(supp)
    supp_edges = CompScienceMeshes.skeleton_fast(supp, 1)
    supp_nodes = CompScienceMeshes.skeleton_fast(supp, 0)

    dir_compl = submesh(!in(dirichlet), bnd_supp)
    dir_compl_edges = CompScienceMeshes.skeleton_fast(dir_compl, 1)
    dir_compl_nodes = CompScienceMeshes.skeleton_fast(dir_compl, 0)

    int_edges = submesh(!in(dir_compl_edges), supp_edges)
    int_nodes = submesh(!in(dir_compl_nodes), supp_nodes)

    Nd_prt = BEAST.nedelec(supp, port_edges)
    Nd_int = BEAST.nedelec(supp, int_edges)
    
    curl_Nd_prt = curl(Nd_prt)
    curl_Nd_int = curl(Nd_int)
    A = assemble(Id, curl_Nd_int, curl_Nd_int, threading=Threading{:single})
    a = -assemble(Id, curl_Nd_int, curl_Nd_prt, threading=Threading{:single}) * x_prt
    
    if length(int_nodes) > 0
        Lg_int = BEAST.lagrangec0d1(supp, int_nodes)
        grad_Lg_int = gradient(Lg_int)
        B = assemble(Id, grad_Lg_int, Nd_int, threading=Threading{:single})
        b = -assemble(Id, grad_Lg_int, Nd_prt, threading=Threading{:single}) * x_prt
        Z = zeros(eltype(b), length(b), length(b))
        S = [A B'; B Z]
        u = S \ [a;b]
    else
        u = A \ a
    end

    x_int = u[1:numfunctions(Nd_int)]
    return x_int, int_edges, Nd_int
end


function extend_2_form(supp, dirichlet, x_prt, port_faces)

    Id = BEAST.Identity()

    bnd_supp = boundary(supp)
    supp_faces = CompScienceMeshes.skeleton_fast(supp, 2)
    supp_edges = CompScienceMeshes.skeleton_fast(supp, 1)
    supp_nodes = CompScienceMeshes.skeleton_fast(supp, 0)

    dir_compl = submesh(!in(dirichlet), bnd_supp)
    dir_compl_faces = CompScienceMeshes.skeleton_fast(dir_compl, 2)
    dir_compl_edges = CompScienceMeshes.skeleton_fast(dir_compl, 1)
    dir_compl_nodes = CompScienceMeshes.skeleton_fast(dir_compl, 0)

    int_faces = submesh(!in(dir_compl_faces), supp_faces)
    int_edges = submesh(!in(dir_compl_edges), supp_edges)
    int_nodes = submesh(!in(dir_compl_nodes), supp_nodes)

    if length(int_nodes) > 0
        @assert length(int_nodes) == 1
        int_edges = int_edges[collect(1:length(int_edges)-length(int_nodes))]
    end

    RT_prt = BEAST.raviartthomas(supp, port_faces)
    RT_int = BEAST.raviartthomas(supp, int_faces)
    Nd_int = BEAST.nedelec(supp, int_edges)

    div_RT_prt = divergence(RT_prt)
    div_RT_int = divergence(RT_int)
    curl_Nd_int = curl(Nd_int)

    A = assemble(Id, div_RT_int, div_RT_int)
    B = assemble(Id, curl_Nd_int, RT_int)

    a = -assemble(Id, div_RT_int, div_RT_prt) * x_prt
    b = -assemble(Id, curl_Nd_int, RT_prt) * x_prt

    Z = zeros(eltype(b), length(b), length(b))
    u = [A B'; B Z] \ [a;b]
    x_int = u[1:end-length(b)]

    return x_int, int_faces, RT_int
end

function dual1forms(Tetrs, Faces, Dir)
    tetrs, bnd, dir, v2t, v2n = dualforms_init(Tetrs, Dir)
    dual1forms_body(Faces, tetrs, bnd, dir, v2t, v2n)
end

function dual2forms(Tetrs, Faces, Dir)
    tetrs, bnd, dir, v2t, v2n = dualforms_init(Tetrs, Dir)
    dual2forms_body(Faces, tetrs, bnd, dir, v2t, v2n)
end



function dual1forms_body(Faces, tetrs, bnd, dir, v2t, v2n)

    @assert dimension(Faces) == 2

    T = coordtype(tetrs)
    bfs = Vector{Vector{Shape{T}}}(undef, length(Faces))
    pos = Vector{vertextype(Faces)}(undef, length(Faces))

    Cells = cells(Faces)
    num_threads = Threads.nthreads()
    Threads.@threads for F in 1:length(Faces)
        Face = Cells[F]
        Chart = chart(Faces, F)

        myid = Threads.threadid()
        myid == 1 && F % 20 == 0 &&
            println("Constructing dual 1-forms: $(F*num_threads) out of $(length(Faces)).")

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
        total_vol = sum(volume(chart(port_edges, edge)) for edge in port_edges)
        nrm = normal(Chart)
        for (i,edge) in enumerate(port_edges)
            cht = chart(port_edges, edge)
            tgt = tangents(center(cht),1)
            vol = volume(cht)
            sgn = sign(dot(nrm, tgt))
            x0[i] = sgn * vol / total_vol
        end

        x23, supp23_int_edges, _ = extend_1_form(supp23, dir23_edges, x0, port_edges)
        x31, supp31_int_edges, _ = extend_1_form(supp31, dir31_edges, x0, port_edges)
        x12, supp12_int_edges, _ = extend_1_form(supp12, dir12_edges, x0, port_edges)

        port1_edges = CompScienceMeshes.union(port_edges, supp31_int_edges)
        port1_edges = CompScienceMeshes.union(port1_edges, supp12_int_edges)
        port2_edges = CompScienceMeshes.union(port_edges, supp12_int_edges)
        port2_edges = CompScienceMeshes.union(port2_edges, supp23_int_edges)
        port3_edges = CompScienceMeshes.union(port_edges, supp23_int_edges)
        port3_edges = CompScienceMeshes.union(port3_edges, supp31_int_edges)

        x1_prt = [x0; x31; x12]
        x2_prt = [x0; x12; x23]
        x3_prt = [x0; x23; x31]

        Nd1_prt = BEAST.nedelec(supp1, port1_edges)
        Nd2_prt = BEAST.nedelec(supp2, port2_edges)
        Nd3_prt = BEAST.nedelec(supp3, port3_edges)

        x1_int, _, Nd1_int = extend_1_form(supp1, dir1_faces, x1_prt, port1_edges)
        x2_int, _, Nd2_int = extend_1_form(supp2, dir2_faces, x2_prt, port2_edges)
        x3_int, _, Nd3_int = extend_1_form(supp3, dir3_faces, x3_prt, port3_edges)

        # inject in the global space
        fn = BEAST.Shape{Float64}[]
        addf!(fn, x1_prt, Nd1_prt, idcs1)
        addf!(fn, x1_int, Nd1_int, idcs1)

        addf!(fn, x2_prt, Nd2_prt, idcs2)
        addf!(fn, x2_int, Nd2_int, idcs2)

        addf!(fn, x3_prt, Nd3_prt, idcs3)
        addf!(fn, x3_int, Nd3_int, idcs3)

        pos[F] = cartesian(CompScienceMeshes.center(Chart))
        bfs[F] = fn
    end

    NDLCCBasis(tetrs, bfs, pos)
end
