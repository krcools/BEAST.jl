using Test

using BEAST
using CompScienceMeshes

V = [
    point(0,0,0),
    point(1,0,0),
    point(0,1,0)]
F = [index(1,2,3)]

m1 = Mesh(V,F)
m2 = CompScienceMeshes.rotate(m1, 0.5pi*point(1,0,0))
m3 = CompScienceMeshes.rotate(m1, 1.0pi*point(1,0,0))

m = weld(m1,m2,m3)

@test numcells(skeleton(m,0)) == 5
@test numcells(m) == 3

Nd = BEAST.nedelec(m)
@test numfunctions(Nd) == 7

function interior(mesh::Mesh, edges=skeleton(mesh,1))
    @assert dimension(mesh) == 2
    @assert vertices(mesh) === vertices(edges)

    C = connectivity(edges, mesh)
    @assert size(C) == (numcells(mesh), numcells(edges))

    nn = vec(sum(abs.(C), dims=1))
    T = CompScienceMeshes.celltype(edges)
    interior_edges = Vector{T}()
    for (i,edge) in pairs(cells(edges))
        nn[i] > 1 && push!(interior_edges, edge)
    end
    Mesh(vertices(mesh), interior_edges)
end

interior_edges = interior(m)
@test numcells(int_edges) == 1
Nd = BEAST.nedelec(m, interior_edges)
@test numfunctions(Nd) == 1
