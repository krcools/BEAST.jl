# build continuous Lagrange elements on a 2D mesh
function lagrangec0(mesh::CompScienceMeshes.AbstractMesh{U,3},
                    nodes::CompScienceMeshes.AbstractMesh{U,1},
                    edges::CompScienceMeshes.AbstractMesh{U,2}; order=1) where {U}  

    T = coordtype(mesh)
    atol = sqrt(eps(T))

    P = vertextype(mesh)
    S = Shape{T}
    F = Vector{S}

    # nodes = skeleton(mesh, 0)
    # edges = skeleton(mesh, 1)

    C02 = connectivity(mesh, nodes, abs); R02 = rowvals(C02); V02 = nonzeros(C02)
    C12 = connectivity(mesh, edges, abs); R12 = rowvals(C12); V12 = nonzeros(C12)

    ne = order-1
    nf = div((order-2)*(order-1), 2)

    nV = length(nodes)
    nE = length(edges) * ne
    nF = length(mesh) * nf
    N = nV + nE + nF

    localspace = LagrangeRefSpace{T,order,3,binomial(2+order,2)}()
    dom = domain(chart(mesh, first(mesh)))
    localdim = numfunctions(localspace, dom)
    
    d = order
    fns = [S[] for n in 1:(nV+nE+nF)]
    zp = zero(SVector{U,T})
    pos = fill(zp, nV+nE+nF)
    for cell in mesh
        cell_ch = chart(mesh, cell)
        V = R02[nzrange(C02,cell)]
        I = V02[nzrange(C02,cell)]
        for (i,v) in zip(I, V)
            vertex_ch = chart(nodes, v)
            gids = [v]
            lids = localindices(_LagrangeGlobalNodesDoFs(d), cell_ch, localspace, i)
            v = globaldofs(vertex_ch, cell_ch, localspace, _LagrangeGlobalNodesDoFs(d))
            α = v[lids, :]
            β = inv(α')
            for i in axes(β,1)
                for j in axes(β,2)
                    isapprox(β[i,j], 0; atol) && continue
                    push!(fns[gids[i]], S(cell, lids[j], β[i,j]))
                end
            end
        end

        order < 2 && continue
        E = R12[nzrange(C12,cell)]
        I = V12[nzrange(C12,cell)]
        for (i,e) in zip(I, E)
            edge_ch = chart(edges, e)
            gids = nV + (e-1)*ne + 1: nV + e*ne
            lids = localindices(_LagrangeGlobalEdgeDoFs(d), cell_ch, localspace, i)
            @assert length(lids) == length(gids)
            v = globaldofs(edge_ch, cell_ch, localspace, _LagrangeGlobalEdgeDoFs(d))
            α = v[lids, :]
            β = inv(α')
            for i in axes(β,1)
                for j in axes(β,2)
                    isapprox(β[i,j], 0; atol) && continue
                    push!(fns[gids[i]], S(cell, lids[j], β[i,j]))
                end
            end
        end

        order < 3 && continue
        F = [cell]
        for (q,f) in enumerate(F)
            face_ch = chart(mesh, f)
            gids = nV + nE + (f-1)*nf + 1 : nV + nE + f*nf
            lids = localindices(_LagrangeGlobalFaceDoFs(d), cell_ch, localspace, 1)
            v = globaldofs(face_ch, cell_ch, localspace, _LagrangeGlobalFaceDoFs(d))
            α = v[lids, :]
            # @show α
            β = inv(α')
            for i in axes(β,1)
                for j in axes(β,2)
                    isapprox(β[i,j], 0; atol) && continue
                    push!(fns[gids[i]], S(cell, lids[j], β[i,j]))
                end
            end
        end 
    end

    for v in nodes pos[v] = cartesian(center(chart(nodes, v))) end
    for e in edges pos[nV + (e-1)*ne + 1: nV + e*ne] .= Ref(cartesian(center(chart(edges, e)))) end
    for f in mesh pos[nV + nE + (f-1)*nf + 1: nV + nE + f*nf] .= Ref(cartesian(center(chart(mesh, f)))) end

    return LagrangeBasis{order,0,localdim}(mesh, fns, pos)
end