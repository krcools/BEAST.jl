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
function buffachristiansen(Γ, γ=mesh(coordtype(Γ),1,3))

    @assert CompScienceMeshes.isoriented(Γ)

    T = coordtype(Γ)

    edges = skeleton(Γ, 1)
    fine = barycentric_refinement(Γ)

    in_interior = interior_tpredicate(Γ)
    on_junction = overlap_gpredicate(γ)

    # first pass to determine the number of functions
    numfuncs = 0
    for edge in cells(edges)
        ch = chart(edges, edge)
        !in_interior(edge) && !on_junction(ch) && continue
        numfuncs += 1
    end

    vtof, vton = vertextocellmap(fine)
    jct_pred = overlap_gpredicate(γ)
    bcs, k = Vector{Vector{Shape{T}}}(undef,numfuncs), 1
    for (i,edge) in enumerate(cells(edges))

        ch = chart(edges, edge)
        !in_interior(edge) && !on_junction(ch) && continue

        # index of edge center in fine's vertexbuffer
        p = numvertices(edges) + i

        v = edge[1]
        n = vton[v]
        F = vec(vtof[v,1:n])
        supp = cells(fine)[F]
        bc1 = buildhalfbc(fine, F, v, p, jct_pred)

        v = edge[2]
        n = vton[v]
        F = vec(vtof[v,1:n])
        supp = cells(fine)[F]
        bc2 = buildhalfbc(fine, F, v, p, jct_pred)
        bc2 = Shape{T}[Shape(s.cellid, s.refid, -s.coeff) for s in bc2]

        bcs[k] = [bc1; bc2]
        k += 1
    end

    return RTBasis(fine, bcs)
end


"""
    buildhalfbc(fine, supp::Array{SVector{3,Int},1}, v, p)
"""
function buildhalfbc(fine, S, v, p, onjunction)

    T = coordtype(fine)

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

    num_ports = length(port_faces)
    port_fluxes = ones(T,num_ports) / num_ports
    @assert length(port_faces) < 3

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
    charges = fill(total_charge/n, n)

    # add the port contributions
    for (f, e, w) in zip(port_faces, port_edges, port_fluxes)
        add!(bf, S[f], e, w)
        charges[f] -= w
    end

    bnd_fluxes = T[]
    if c_on_boundary
        if num_junctions == 0
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

    return bf
end
