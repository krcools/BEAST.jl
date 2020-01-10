using SparseArrays

function edge_values(space, m)
    rs = refspace(space)
    supp = geometry(space)
    edgs = skeleton(supp,1)
    Conn = copy(transpose(connectivity(edgs, supp, identity)))
    Vals = zeros(scalartype(space), length(supp), length(edgs))
    for (i,sh) in enumerate(space.fns[m])
        tet_id = sh.cellid
        tet = chart(supp, cells(supp)[tet_id])
        for k in nzrange(Conn, tet_id)
            edg_id = rowvals(Conn)[k]
            # loc_id = abs(nonzeros(Conn)[k])
            # edge = CompScienceMeshes.edges(tet)[loc_id]
            edge = chart(edgs, cells(edgs)[edg_id])
            ctr = center(edge)
            tgt = tangents(ctr,1)
            ctr1 = neighborhood(tet, carttobary(tet, cartesian(ctr)))
            vals = rs(ctr1)
            # @show vals
            val = vals[sh.refid].value
            Vals[tet_id, edg_id] += dot(sh.coeff * val, tgt)
        end
    end
    return Vals
end

function check_edge_values(EV)
    for i in axes(EV,2)
        ev = [x for x in EV[:,i] if abs(x) > 1e-8]
        uev = unique(x -> round(x, digits=5), ev)
        if !(length(uev) <= 1)
            @show i
            @show uev
            error("stop")
        end
    end
end
