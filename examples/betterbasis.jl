using CompScienceMeshes, BEAST
Base.getindex(m::CompScienceMeshes.AbstractMesh, i::Int) = cells(m)[i]
cellvertices(m::CompScienceMeshes.AbstractMesh, i::Int) = vertices(m)[cells(m)[i]]

m = meshsphere(1.0, 0.35)
# m = meshcuboid(1.0, 1.0, 0.2, 0.2)

X = raviartthomas(m)
Y = buffachristiansen(m, ibscaled=true)

Id = BEAST.Identity()
Nx = BEAST.NCross()

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=1.0)
G = -assemble(Id, Y, Y)

b = assemble((n×E)×n, Y)
x = G \ b

# x = zeros(Float64,numfunctions(Y))
# x[20] = 1
fcr, geo = facecurrents(x, Y)

import PlotlyJS
using LinearAlgebra
PlotlyJS.plot(patch(geo.mesh, norm.(fcr)))
# PlotlyJS.plot(patch(geo.mesh, rand(numcells(geo))))

edges = skeleton(m,1)
verts = skeleton(m,0)

vertex_list = [v[1] for v in cells(verts)]
Q = lagrangec0d1(m, vertex_list, Val{3})
Z = curl(Q)

Gzy = assemble(Id, Z, Y)

# vtoc, nc = vertextocellmap(edges)
D = connectivity(verts, edges)
@assert size(D) == (numcells(edges), numcells(verts))

los = geometry(Y)
# Choose a vertex:
for n_vert in 1:numcells(verts)
    @assert cellvertices(verts, n_vert)[1] ≈ Q.pos[n_vert]
# find the connected edges
    # n_edges = vtoc[n_vert,1:nc[n_vert]]
    using SparseArrays
    rows = rowvals(D)
    vals = nonzeros(D)
    for k in nzrange(D,n_vert)
        n_edge = rows[k]
        cht = chart(edges, edges[n_edge])
        ctr = cartesian(center(cht))
        vt1 = cht.vertices[1]
        lvt = los.mesh.vertices[numvertices(m)+n_edge]
        pos = Y.pos[n_edge]
        @assert norm(lvt - pos) ≤ 1e-8
        @assert norm(cross(pos-vt1, ctr-vt1)) ≤ 1e-8
    end
end

n_vert = 13
x = [float(i == n_vert) for i in 1:numfunctions(Q)]
fcr1, geo1 = facecurrents(x, Q)
colors = getindex.(fcr1,1)
PlotlyJS.plot(patch(m, colors))

divY = divergence(Y)
fcr2, geo2 = facecurrents(rand(numfunctions(divY)), divY)
colors = getindex.(fcr2,1)
PlotlyJS.plot(patch(geo2.mesh, colors))

for n_fn in 1:numfunctions(divY)
    n_fn = 21
    charges = zeros(numcells(los))
    for shape in divY.fns[n_fn]
        f = shape.cellid
        ch = chart(los, los[f])
        charges[f] += shape.coeff #* volume(ch)
    end
    nzs = (charges[(charges .≈ 0) .== false])
    xtr = extrema(nzs)
    @assert all( Vector{Bool}(nzs .≈ xtr[1]) .| Vector{Bool}(nzs .≈ xtr[2]))
end
