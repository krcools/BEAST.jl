export buffachristiansen

"""
Move the s-th element right after the d-th
"""
function move_after!(p,n,s,d)

    n[d] == s && return

    t1 = n[d]
    t2 = n[s]
    t3 = p[s]

    n[d] = s
    n[s] = t1
    p[s] == 0 || (n[p[s]] = t2)

    p[s] = d
    t1 == 0 || (p[t1] = s)
    t2 == 0 || (p[t2] = t3)
end

function move_before!(p, n, s, d)
    p[d] == s && return

    t1 = p[d]
    t2 = p[s]
    t3 = n[s]

    p[d] = s
    p[s] = t1
    n[s] == 0 || (p[n[s]] = t2)

    n[s] = d
    t1 == 0 || (n[t1] = s)
    t2 == 0 || (n[t2] = t3)
end

function sortneighbors(a, pred)

    n = collect(2:length(a)+1); n[end] = 0
    p = collect(0:length(a)-1); p[1] = 0

    last = 1
    while last != 0
        cand = n[last]
        while cand != 0
            pred(a[last], a[cand]) && break
            cand = n[cand]
        end
        cand == 0 && break
        move_after!(p, n, cand, last)
        last = cand
    end

    first = 1
    while last != 0
        cand = n[last]
        while cand != 0
            pred(a[cand], a[first]) && break
            cand = n[cand]
        end
        cand == 0 && break
        move_before!(p, n, cand, first)
        first = cand
    end

    b = similar(a)
    i, j = last, length(n)
    while true
        b[j] = a[i]
        i = p[i]
        i == 0 && break
        j -= 1
    end

    return b

end

isclosed(a, pred) = length(a)>2 && pred(a[end], a[1])


"""
    buffachristiansen(Γ, γ)

Construct the set of Buffa-Christiansen functions subject to mesh Γ and only
enforcing zero normal components on ∂Γ ∖ γ.
"""
function buffachristiansen(Γ, γ=mesh(coordtype(Γ),1,3); ibscaled=false, sort=:spacefillingcurve, edges=:all)

    @assert CompScienceMeshes.isoriented(Γ)

    T = coordtype(Γ)
    P = vertextype(Γ)

    in_interior = interior_tpredicate(Γ)
    on_junction = overlap_gpredicate(γ)

    # edges = skeleton(Γ, 1; sort)

    if edges == :all
        edges = skeleton(Γ, 1; sort)
        edges = submesh(edges) do edges, edge
            ch = chart(edges, edge)
            !in_interior(edges, edge) && !on_junction(ch) && return false
            return true
        end
    end
    @assert edges isa CompScienceMeshes.AbstractMesh
    @assert dimension(edges) == 1

    fine = if ibscaled
        CompScienceMeshes.lineofsight_refinement(Γ)
    else
        barycentric_refinement(Γ; sort)
    end


    # # first pass to determine the number of functions
    # numfuncs = 0
    # for edge in edges
    #     ch = chart(edges, edge)
    #     !in_interior(edge) && !on_junction(ch) && continue
    #     numfuncs += 1
    # end

    vtof, vton = vertextocellmap(fine)
    jct_pred = overlap_gpredicate(γ)
    bcs, k = Vector{Vector{Shape{T}}}(undef,length(edges)), 1
    pos = Vector{P}(undef,length(edges))
    for (i,e) in enumerate(edges)

        edge = CompScienceMeshes.indices(edges,e)
        ch = chart(edges, e)
        ln = volume(ch)
        # !in_interior(edge) && !on_junction(ch) && continue

        # index of edge center in fine's vertexbuffer
        v  = edge[1]
        n = vton[v]
        F = vtof[v,1:n]
        supp = fine[F]
        bnd = boundary(supp)
        bnd_nodes = skeleton(bnd,0)
        p = 0
        for bn in bnd_nodes
            node = CompScienceMeshes.indices(bnd_nodes, bn)
            vert1 = vertices(Γ)[edge[1]]
            vert2 = vertices(Γ)[edge[2]]
            vert = vertices(fine)[node[1]]
            dist = norm(vert-vert1)
            node[1] == edge[1] && continue
            if dot(vert-vert1,vert2-vert1) ≈ dist*ln
                p = node[1]
                break
            end
        end
        @assert p != 0

        # p = numvertices(edges) + i
        pos[k] = vertices(fine)[p]

        # sanity check
        edge_ctr = cartesian(center(ch))
        fine_vtx = vertices(fine)[p]
        @assert all(carttobary(ch, edge_ctr) .≥ 0)
        @assert norm(cross(edge_ctr - ch.vertices[1], fine_vtx - ch.vertices[1])) ≤ 1e-8

        v = edge[1]
        n = vton[v]
        F = vec(vtof[v,1:n])
        supp = cells(fine)[F]
        bc1 = buildhalfbc(fine, F, v, p, jct_pred, ibscaled)

        v = edge[2]
        n = vton[v]
        F = vec(vtof[v,1:n])
        supp = cells(fine)[F]
        bc2 = buildhalfbc(fine, F, v, p, jct_pred, ibscaled)
        bc2 = Shape{T}[Shape(s.cellid, s.refid, -s.coeff) for s in bc2]

        bcs[k] = [bc1; bc2]
        k += 1
    end

    return RTBasis(fine, bcs, pos)
end


"""
    buildhalfbc(fine, supp::Array{SVector{3,Int},1}, v, p)
"""
function buildhalfbc(fine, S, v, p, onjunction, ibscaled)

    T = coordtype(fine)

    @assert v != 0
    @assert p != 0
    @assert length(S) >= 2
    @assert mod(length(S), 2) == 0
    @assert all(0 .< S .<= numcells(fine))
    @assert all([v in cells(fine)[s] for s in S])

    bf = Shape{T}[]

    n = length(S)
    N = n ÷ 2
    modn(i) = mod1(i,n)

    # sort the fine faces making up the domain
    share_edge(i,j) = length(intersect(cells(fine)[i],cells(fine)[j])) == 2
    S = sortneighbors(S, share_edge)
    @assert all(0 .< S .<= numcells(fine))

    port_faces = findall(i -> p in cells(fine)[i], S)
    port_edges = [something(findfirst(isequal(v), cells(fine)[S[i]]),0) for i in port_faces]

    c_on_boundary = !share_edge(S[end], S[1]) || length(S) == 2
    e_on_boundary = length(port_faces) == 1
    @assert c_on_boundary || !e_on_boundary

    # if there are two ports (c must be interior), we want to iterate
    # the support starting with one and ending with the other. This
    # might require reordering the ports we just discovered.
    if !e_on_boundary
        @assert length(port_faces) == 2 "$port_faces, $n"
        if port_faces[2] == modn(port_faces[1]+1)
            reverse!(port_faces)
            reverse!(port_edges)
        end
        @assert port_faces[1] == modn(port_faces[2]+1)
    end

    port_lengths = [begin
            face = cells(fine)[S[j]]
            v1 = vertices(fine)[face[mod1(i+1,3)]]
            v2 = vertices(fine)[face[mod1(i+2,3)]]
            norm(v2-v1)
        end for (j,i) in zip(port_faces, port_edges)]
    num_ports = length(port_faces)
    port_fluxes = if ibscaled
        total_length = sum(port_lengths)
        port_lengths / total_length
    else
        ones(T,num_ports) / num_ports
    end
    @assert sum(port_fluxes) ≈ 1
    @assert length(port_faces) < 3

    # Detect the boundary edges
    bnd_faces = Int[]
    bnd_edges = Int[]
    if c_on_boundary
        bnd_faces = [1,n]
        p1 = x -> x in cells(fine)[S[2]] && x != v
        p2 = x -> x in cells(fine)[S[n-1]] && x != v
        e1 = findfirst(p1, cells(fine)[S[1]])
        e2 = findfirst(p2, cells(fine)[S[n]])
        bnd_edges = [e1, e2]
    end

    # Detect the junction edges
    num_junctions = 0
    jct_idxes = Int[]
    if c_on_boundary
        for i in 1:length(bnd_faces)
            f, e = bnd_faces[i], bnd_edges[i]
            a = vertices(fine)[cells(fine)[S[f]][mod1(e+1,3)]]
            b = vertices(fine)[cells(fine)[S[f]][mod1(e+2,3)]]
            seg = simplex(a,b)
            if onjunction(seg)
                push!(jct_idxes, i)
                #push!(jct_edges, e)
            end
        end
        num_junctions = length(jct_idxes)
    end

    # This charge needs to be compensated by interior divergence
    total_charge = (!c_on_boundary || num_junctions == 2) ? 1 : 0
    charges = if ibscaled
        face_areas = [volume(chart(fine, s)) for s in S]
        face_areas / sum(face_areas)
    else
        fill(total_charge/n, n)
    end
    @assert sum(charges) ≈ 1 || isapprox(sum(charges), 0, atol=1e-8)

    # add the port contributions
    for (f, e, w) in zip(port_faces, port_edges, port_fluxes)
        add!(bf, S[f], e, w)
        charges[f] -= w
    end

    bnd_fluxes = T[]
    if c_on_boundary
        if num_junctions == 0
            ibscaled && error("IB scaled BCs for mesh with boundary not implemented!")
            bnd_fluxes = T[-(n-port_faces[2])/n, -port_faces[2]/n]
        elseif num_junctions == 1
            bnd_fluxes = -ones(T,2)
            bnd_fluxes[jct_idxes[1]] = 0
        else
            bnd_fluxes = zeros(T,2)
        end
    end

    @assert length(bnd_faces) == length(bnd_edges) == length(bnd_fluxes)

    for (f, e, w) in zip(bnd_faces, bnd_edges, bnd_fluxes)
        add!(bf, S[f], e, w)
        charges[f] -= w
    end

    # The remaining fluxes can be computed based on port and boundary
    # influx, and the knowledge of the charge density by recursion.
    # Iterate beginning with port 0:
    start_face = c_on_boundary ? 1 : port_faces[1]
    for i in 1 : n-1
        # cyclic port_face[1]-based indexing
        j0 = modn(i + start_face - 1)
        j1 = modn(j0+1)
        # get the indices of the faces in fine.faces
        f0, f1 = S[j0], S[j1]
        # what are the local indices of the common edge?
        e0, e1 = getcommonedge(cells(fine)[f0], cells(fine)[f1])
        # add the two half Raviart-Thomas shape functions
        add!(bf, S[j0], abs(e0), +charges[j0])
        add!(bf, S[j1], abs(e1), -charges[j0])
        # update the charge bookkeeping
        charges[j1] += charges[j0]
    end

    # Add a multiple of the zero-divergence pattern such that the result
    # is orthogonal to the zero-divergence pattern
    if ibscaled
        β = zero(T)
        for shape in bf
            f = shape.cellid
            @assert f in S
            face = cells(fine)[f]
            i = something(findfirst(==(v),face), 0)
            @assert i != 0
            ch = chart(fine, f)
            area = volume(ch)
            qps = quadpoints(ch, 3)
            @assert sum(w for (p,w) in qps) ≈ area
            for (p,w) in qps
                vals = RTRefSpace{T}()(p)
                bfp = shape.coeff * vals[shape.refid].value
                i1 = mod1(i+1,3)
                i2 = mod1(i+2,3)
                dfp = 0.5/area * (ch.vertices[i2] - ch.vertices[i1])
                @assert dfp ≈ vals[i1].value - vals[i2].value
                β += w * dot(bfp, dfp)
            end
        end

        γ = zero(T)
        for f in S
            face = cells(fine)[f]
            i = something(findfirst(==(v), face), 0)
            @assert i != 0
            ch = chart(fine, f)
            area = volume(ch)
            vct = ch.vertices[mod1(i+2,3)] - ch.vertices[mod1(i+1,3)]
            γ += 0.25/area * dot(vct,vct)
        end

        α = -β/γ
        for f in S
            face = cells(fine)[f]
            i = something(findfirst(==(v), face), 0)
            @assert i != 0
            add!(bf, f, mod1(i+1,3), α)
            add!(bf, f, mod1(i+2,3), -α)
        end

        β = zero(T)
        for shape in bf
            f = shape.cellid
            @assert f in S
            face = cells(fine)[f]
            i = something(findfirst(==(v),face), 0)
            @assert i != 0
            ch = chart(fine, f)
            area = volume(ch)
            qps = quadpoints(ch, 3)
            @assert sum(w for (p,w) in qps) ≈ area
            for (p,w) in qps
                vals = RTRefSpace{T}()(p)
                bfp = shape.coeff * vals[shape.refid].value
                i1 = mod1(i+1,3)
                i2 = mod1(i+2,3)
                dfp = 0.5/area * (ch.vertices[i2] - ch.vertices[i1])
                @assert dfp ≈ vals[i1].value - vals[i2].value
                β += w * dot(bfp, dfp)
            end
        end
        @assert isapprox(β,0, atol=1e-8)
    end

    return bf
end


using LinearAlgebra
function buildhalfbc2(patch, port, dirichlet, prt_fluxes)

    # println()

    edges = skeleton(patch,1)
    verts = skeleton(patch,0)
    bndry = boundary(patch)

    in_bndry = in(bndry)
    dirbnd = submesh(dirichlet) do m,edge
        in_bndry(m,edge)
        # edge in cells(bndry) && return true
        # reverse(edge) in cells(bndry) && return true
        # return false
    end
    # @show numcells(patch)
    # @show numcells(dirbnd)
    @assert numcells(dirbnd) ≤ 4

    bnd_dirbnd = boundary(dirbnd)
    nodes_dirbnd = skeleton(dirbnd,0)
    in_bnd_dirbnd = in(bnd_dirbnd)
    int_nodes_dirbnd = submesh(nodes_dirbnd) do m,node
        # node in cells(bnd_dirbnd) && return false
        # return true
        return !in_bnd_dirbnd(m, node)
    end
    # @show numcells(int_nodes_dirbnd)
    @assert numcells(int_nodes_dirbnd) ≤ 2

    int_pred = interior_tpredicate(patch)
    num_edges_on_port = 0
    num_edges_on_dirc = 0
    # @show numcells(edges)

    # for edge in cells(edges)
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

    in_port = in(port)
    in_dirbnd = in(dirbnd)
    int_edges = submesh(edges) do m,edge

        in_port(m,edge) && (num_edges_on_port+=1; return false)
        in_dirbnd(m,edge) && (num_edges_on_dirc+=1; return true)
        !int_pred(m,edge) && return false
        return true

        # (edge in cells(port)) && (num_edges_on_port+=1 ; return false)
        # (reverse(edge) in cells(port)) && (num_edges_on_port+=1 ; return false)
        # edge in cells(dirbnd) && (num_edges_on_dirc+=1; return true)
        # reverse(edge) in cells(dirbnd) && (num_edges_on_dirc+=1; return true)
        # (!int_pred(edge)) && return false
        # return true
    end
    # println()
    # for edge in cells(int_edges)
    #     println(edge)
    # end
    # @show numcells(int_edges)
    # @show num_edges_on_port
    # @show num_edges_on_dirc

    bnd_verts = skeleton(bndry,0)
    prt_verts = skeleton(port,0)
    # function isonboundary(node)
    #     node in cells(bnd_verts)
    # end
    # int_verts = submesh(!isonboundary, verts)
    # dirichlet_nodes = skeleton(dirichlet,0)
    in_int_nodes_dirbnd = in(int_nodes_dirbnd)
    in_bnd_verts = in(bnd_verts)
    in_prt_verts = in(prt_verts)
    int_verts = submesh(verts) do m,node
        in_int_nodes_dirbnd(m,node) && return true
        in_bnd_verts(m,node) && return false
        in_prt_verts(m,node) && return false
        return true
        # node in cells(int_nodes_dirbnd) && return true
        # node in cells(bnd_verts) && return false
        # node in cells(prt_verts) && return false
        # return true
    end
    # @show numcells(int_verts)

    RT_int = raviartthomas(patch, cellpairs(patch, int_edges))
    RT_prt = raviartthomas(patch, cellpairs(patch, port))
    vertex_list = [c[1] for c in cells(int_verts)]
    L0_int = lagrangec0d1(patch, vertex_list, Val{3})


    Id = BEAST.Identity()
    D = assemble(Id, divergence(RT_int), divergence(RT_int))
    Q = assemble(Id, divergence(RT_int), divergence(RT_prt))
    d = -Q * prt_fluxes

    # @show numfunctions(L0_int)
    @assert numfunctions(L0_int) in [1,2]
    C = assemble(Id, curl(L0_int), RT_int)
    curl_L0_int = curl(L0_int)
    c = real(assemble(Id, curl_L0_int, RT_prt)) * prt_fluxes

    mD = Matrix(D)
    x1 = pinv(mD) * d
    N = nullspace(mD)
    # @show size(N)
    # @show size(C)
    p = (C*N) \ (c - C*x1)
    x = x1 + N*p

    @assert D*x ≈ d atol=1e-8
    if !isapprox(C*x, c, atol=1e-8)
        @show norm(C*x-c)
        error("error")
    end

    return RT_int, RT_prt, x, prt_fluxes

end

# Use the algebraic construction also used in dual3d
function buffachristiansen2(Faces::CompScienceMeshes.AbstractMesh{U,D1,T}) where {U,D1,T}

    faces = barycentric_refinement(Faces)
    Edges = skeleton(Faces,1)
    Bndry = boundary(Faces)
    inBndry = in(Bndry)
    Edges = submesh(Edges) do m,Edge
        return !inBndry(m,Edge)
        # Edge in cells(Bndry) && return false
        # reverse(Edge) in cells(Bndry) && return false
        # return true
    end

    #T = Float64
    bfs = Vector{Vector{Shape{T}}}(undef, numcells(Edges))
    pos = Vector{vertextype(Faces)}(undef, numcells(Edges))
    dirichlet = boundary(faces)
    for (E,Edge) in enumerate(Edges)
        EdgeInds = CompScienceMeshes.indices(Edges, Edge)
        # @show Edge
        bfs[E] = Vector{Shape{T}}()

        pos[E] = cartesian(center(chart(Edges,Edge)))
        # port_vertex_idx = numvertices(Faces) + E
        port_vertex_idx = argmin(norm.(vertices(faces) .- Ref(pos[E])))
        # pos[E] = vertices(faces)[port_vertex_idx]
        # @show carttobary(chart(Edges,Edge),pos[E])
        # @show pos[E]
        # @show (Faces.vertices[Edge[1]]+Faces.vertices[Edge[2]])/2

        # Build the plus-patch
        ptch_vert_idx = EdgeInds[1]
        ptch_face_idcs = [i for (i,face) in enumerate(cells(faces)) if ptch_vert_idx in face]
        patch = Mesh(vertices(faces), cells(faces)[ptch_face_idcs])
        patch_bnd = boundary(patch)
        # @show numcells(patch_bnd)
        port = Mesh(vertices(faces), filter(c->port_vertex_idx in c, cells(patch_bnd)))

        # @show numcells(patch)
        # @show numcells(port)

        # @assert numcells(patch) >= 6
        @assert numcells(port) == 2

        prt_fluxes = [0.5, 0.5]
        RT_int, RT_prt, x_int, x_prt = buildhalfbc2(patch, port, dirichlet, prt_fluxes)

        for (m,bf) in enumerate(RT_int.fns)
            for sh in bf
                cellid = ptch_face_idcs[sh.cellid]
                BEAST.add!(bfs[E], cellid, sh.refid, sh.coeff * x_int[m])
            end
        end

        for (m,bf) in enumerate(RT_prt.fns)
            for sh in bf
                cellid = ptch_face_idcs[sh.cellid]
                BEAST.add!(bfs[E],cellid, sh.refid, sh.coeff * x_prt[m])
            end
        end

        # Build the minus-patch
        ptch_vert_idx = EdgeInds[2]
        ptch_face_idcs = [i for (i,face) in enumerate((faces)) if ptch_vert_idx in CompScienceMeshes.indices(faces, face)]
        patch = Mesh(vertices(faces), cells(faces)[ptch_face_idcs])
        port = Mesh(vertices(faces), filter(c->port_vertex_idx in c, cells(boundary(patch))))
        # @assert numcells(patch) >= 6
        @assert numcells(port) == 2

        prt_fluxes = [-0.5, -0.5]
        RT_int, RT_prt, x_int, x_prt = buildhalfbc2(patch, port, dirichlet, prt_fluxes)

        for (m,bf) in enumerate(RT_int.fns)
            for sh in bf
                cellid = ptch_face_idcs[sh.cellid]
                BEAST.add!(bfs[E], cellid, sh.refid, sh.coeff * x_int[m])
            end
        end

        for (m,bf) in enumerate(RT_prt.fns)
            for sh in bf
                cellid = ptch_face_idcs[sh.cellid]
                BEAST.add!(bfs[E],cellid, sh.refid, sh.coeff * x_prt[m])
            end
        end

    end

    RTBasis(faces, bfs, pos)
end

# Extend into both dual faces in a single call to buildhalfbc2 (aka extend_2_form)
function buffachristiansen3(Faces::CompScienceMeshes.AbstractMesh{U,D1,T}) where {U,D1,T}

    faces = barycentric_refinement(Faces)
    Edges = skeleton(Faces,1)
    Bndry = boundary(Faces)
    inBndry = in(Bndry)
    Edges = submesh(Edges) do m,Edge
        inBndry(m,Edge) && return false
        # reverse(Edge) in cells(Bndry) && return false
        return true
    end

    #T = Float64
    bfs = Vector{Vector{Shape{T}}}(undef, numcells(Edges))
    pos = Vector{vertextype(Faces)}(undef, numcells(Edges))
    dirichlet = boundary(faces)
    for (E,Edge) in enumerate(Edges)
        bfs[E] = Vector{Shape{T}}()

        # port_vertex_idx = numvertices(Faces) + E
        # pos[E] = vertices(faces)[port_vertex_idx]
        pos[E] = cartesian(center(chart(Edges,Edge)))
        # port_vertex_idx = numvertices(Faces) + E
        port_vertex_idx = argmin(norm.(vertices(faces) .- Ref(pos[E])))

        # Build the dual support
        EdgeInds = CompScienceMeshes.indices(Edges, Edge)
        ptch_face_idcs1 = [i for (i,face) in enumerate(faces) if EdgeInds[1] in CompScienceMeshes.indices(faces,face)]
        ptch_face_idcs2 = [i for (i,face) in enumerate(faces) if EdgeInds[2] in CompScienceMeshes.indices(faces,face)]
        patch1 = Mesh(vertices(faces), cells(faces)[ptch_face_idcs1])
        patch2 = Mesh(vertices(faces), cells(faces)[ptch_face_idcs2])
        patch_bnd = boundary(patch1)
        port = Mesh(vertices(faces), filter(c->port_vertex_idx in c, cells(patch_bnd)))
        patch = CompScienceMeshes.union(patch1, patch2)
        # @show numcells(patch_bnd)
        # @show numcells(patch)
        # @show numcells(port)
        # @assert numcells(skeleton(patch,1))+2 == numcells(skeleton(patch1,1)) + numcells(skeleton(patch2,1))

        @assert numcells(patch) >= 6
        @assert numcells(port) == 2
        prt_fluxes = [0.5, -0.5]
        RT_int, RT_prt, x_int, x_prt = buildhalfbc2(patch, port, dirichlet, prt_fluxes)

        ptch_face_idcs = vcat(ptch_face_idcs1, ptch_face_idcs2)
        for (m,bf) in enumerate(RT_int.fns)
            for sh in bf
                cellid = ptch_face_idcs[sh.cellid]
                BEAST.add!(bfs[E], cellid, sh.refid, sh.coeff * x_int[m])
            end
        end

        for (m,bf) in enumerate(RT_prt.fns)
            for sh in bf
                cellid = ptch_face_idcs[sh.cellid]
                BEAST.add!(bfs[E],cellid, sh.refid, sh.coeff * x_prt[m])
            end
        end

        # # Build the minus-patch
        # ptch_vert_idx = Edge[2]
        # ptch_face_idcs = [i for (i,face) in enumerate(cells(faces)) if ptch_vert_idx in face]
        # patch = Mesh(vertices(faces), cells(faces)[ptch_face_idcs])
        # port = Mesh(vertices(faces), filter(c->port_vertex_idx in c, cells(boundary(patch))))
        # @assert numcells(patch) >= 6
        # @assert numcells(port) == 2
        # RT_int, RT_prt, x_int, x_prt = buildhalfbc2(patch, port, nothing)
    
        # for (m,bf) in enumerate(RT_int.fns)
        #     for sh in bf
        #         cellid = ptch_face_idcs[sh.cellid]
        #         BEAST.add!(bfs[E], cellid, sh.refid, -1.0 * sh.coeff * x_int[m])
        #     end
        # end
    
        # for (m,bf) in enumerate(RT_prt.fns)
        #     for sh in bf
        #         cellid = ptch_face_idcs[sh.cellid]
        #         BEAST.add!(bfs[E],cellid, sh.refid, -1.0 * sh.coeff * x_prt[m])
        #     end
        # end

    end

    RTBasis(faces, bfs, pos)
end
