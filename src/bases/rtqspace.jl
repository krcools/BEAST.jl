struct RTQSpace{T,M,P} <: Space{T}
    geo::M
    fns::Vector{Vector{Shape{T}}}
    pos::Vector{P}
end

# RTQSpace(g::M, fns::Vector{})

function positions(s::RTQSpace) s.pos end 
function refspace(s::RTQSpace{T}) where {T} RTQuadRefSpace{T}() end
function subset(rt::RTQSpace,I) RTQSpace(rt.geo, rt.fns[I], rt.pos[I]) end

function raviartthomas(
    mesh::CompScienceMeshes.QuadMesh{T},
    edges::CompScienceMeshes.AbstractMesh{3,2,T},
    connectivity, orientations) where {T<:Any}

    fns = Vector{Vector{Shape{T}}}(undef, length(edges))
    pos = Vector{SVector{3,T}}(undef, length(edges))

    rows = rowvals(connectivity)
    vals = nonzeros(connectivity)
    for (i,edge) in enumerate(edges)
        σ = orientations[i]
        fn = map(zip(nzrange(connectivity, i),(σ,-σ))) do (j, α)
            Shape{T}(rows[j], abs(vals[j]), α)
        end
        fns[i] = fn
        pos[i] = cartesian(CompScienceMeshes.center(chart(edges, edge)))
    end

    return RTQSpace(mesh, fns, pos)
end

function raviartthomas(
    mesh::CompScienceMeshes.QuadMesh{T},
    edges::CompScienceMeshes.AbstractMesh{3,2,T},
    orientations::Vector{Bool}) where {T<:Any}

    conn = connectivity(edges, mesh)
    return raviartthomas(mesh, edges, conn, orientations)
end

function raviartthomas(
    mesh::CompScienceMeshes.QuadMesh{T},
    edges::CompScienceMeshes.AbstractMesh{3,2,T}) where {T}

    c = connectivity(edges, mesh, identity)
    o = ones(length(edges))

    vals = nonzeros(c)
    rows = rowvals(c)
    for i in axes(c,2)
        k = first(nzrange(c,i))
        vals[k] < 0 && (o[i] = -1)
    end

    return raviartthomas(mesh, edges, c, o)
end

function raviartthomas(mesh::CompScienceMeshes.QuadMesh{T}) where {T}
    edges = skeleton(mesh,1)
    bnd = boundary(mesh)
    edges_int = submesh(!in(bnd), edges)
    raviartthomas(mesh, edges_int)
end

@testitem "RTQSpace construction" begin
    using CompScienceMeshes

    m = CompScienceMeshes.meshrectangle(2.0, 2.0, 1.0; element=:quadrilateral)
    edges = skeleton(m, 1)
    edges_bnd = boundary(m)
    # @show length(edges_bnd)
    @test length(edges_bnd) == 8
    pred = !in(edges_bnd)
    edges_int = submesh(pred,  edges)

    c = CompScienceMeshes.connectivity(edges_int, m, identity)
    @test size(c) == (length(m), length(edges_int))
    o = ones(length(edges_int))
    s = raviartthomas(m, edges_int, c, o)
    @test numfunctions(s) == 4
end

@testitem "RTQSpace assembly data" begin
    using CompScienceMeshes
    m = CompScienceMeshes.meshrectangle(2.0, 2.0, 1.0; element=:quadrilateral)
    edges = skeleton(m, 1)
    edges_bnd = boundary(m)
    pred = !in(edges_bnd)
    edges_int = submesh(pred,  edges)
    c = CompScienceMeshes.connectivity(edges_int, m, identity)
    o = ones(length(edges_int))
    s = raviartthomas(m, edges_int, c, o)

    num_cells = numcells(m)
    num_bfs = numfunctions(s)
    r = refspace(s)
    dom = domain(chart(m, first(m)))
    num_refs = numfunctions(r, dom)
    celltonum = BEAST.make_celltonum(num_cells, num_refs, num_bfs, s)

    els, ad, a2g = BEAST.assemblydata(s)
end