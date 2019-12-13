using CompScienceMeshes
using BEAST

tetrs = CompScienceMeshes.tetmeshsphere(1.0, 0.20)
@show numcells(tetrs)

bndry = boundary(tetrs)
edges = skeleton(tetrs, 1)

bndry_edges = [sort(c) for c in cells(skeleton(bndry, 1))]
function is_interior(edge)
    !(sort(edge) in bndry_edges)
end

interior_edges = submesh(is_interior, edges)
@assert numcells(interior_edges) + numcells(skeleton(bndry,1)) == numcells(edges)

X = BEAST.nedelecc3d(tetrs, interior_edges)
@assert numfunctions(X) == numcells(interior_edges)

Id = BEAST.Identity()
Y = curl(X)
A1 = assemble(Id, Y, Y)
A2 = assemble(Id, X, X)
A = A1 - A2

using LinearAlgebra
f = BEAST.ScalarTrace(x -> point(1,0,0) * exp(-norm(x)^2/4))
b = assemble(f, X)

u = A \ b

using Plots

pts = [point(0,t,0) for t in range(-2,2,length=200)];
vals = BEAST.grideval(pts, u, X)
plot(getindex.(real(vals),1), m=2)
