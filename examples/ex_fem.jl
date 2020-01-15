using CompScienceMeshes
using BEAST

using SparseArrays
function isdivconforming(space)

    geo = geometry(space)
    mesh = geo

    edges = skeleton(mesh,1)
    D = connectivity(edges, mesh, abs)
    rows = rowvals(D)
    vals = nonzeros(D)

    Flux = zeros(Float64, numcells(mesh),3)
    TotalFlux = zeros(Float64, numcells(edges), numfunctions(space))

    for i in 1 : numfunctions(space)

        fill!(Flux, 0)
        bfs = BEAST.basisfunction(space,i)
        for bf in bfs
            c = bf.cellid
            e = bf.refid
            x = bf.coeff
            Flux[c,e] += x
        end

        for E in 1 : numcells(edges)
            for j in nzrange(D,E)
                F = rows[j]
                e = vals[j]
                @assert 1 <= e <= 3
                TotalFlux[E,i] += Flux[F,e]
            end
        end
    end

    I = findall(abs.(TotalFlux) .> 1e-6)
    @show length(I)
    return TotalFlux, I
end

tetrs = CompScienceMeshes.tetmeshsphere(1.0, 0.15)
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
# f = BEAST.ScalarTrace(x -> point(1,0,0) * exp(-norm(x)^2/4))
f = BEAST.ScalarTrace(x -> point(1,0,0))
b = assemble(f, X)

u = A \ b

using Plots

pts = [point(0,t,0) for t in range(-2,2,length=200)];
vals = BEAST.grideval(pts, u, X)
vals2 = BEAST.grideval(pts, u, Y)
vals2 = [point(0,1,0) × val for val in vals2]
plot(getindex.(real(vals),1), m=2)
plot(norm.(vals2), m=2)

# Inspect the value of a single basis function
p = 120
tet = chart(tetrs, cells(tetrs)[p])
nbd = neighborhood(tet, [0.5, 0.5, 0.0])
edg = simplex(tet.vertices[1], tet.vertices[2])
tgt = normalize(edg.vertices[2] - edg.vertices[1])

vals = refspace(X)(nbd)
@assert dot(tgt,vals[1].value) ≈ 1/volume(edg) atol=1e-8

# check also the value of the raviart-thomas 3d in the center
# of their defining face
Z = BEAST.nedelecd3d(tetrs)
@assert numfunctions(Z) == numcells(skeleton(tetrs,2))

fce = simplex(tet.vertices[2], tet.vertices[3], tet.vertices[4])
nbd_rt = neighborhood(tet, [0.0, 1/3, 1/3])
vals_rt = refspace(Z)(nbd_rt)

# carttobary(fce, cartesian(nbd_rt))

#
dot(vals_rt[1].value, normal(fce)) * volume(fce)
o, x, y, z = euclidianbasis(3)
@assert volume(simplex(x,y,z,o)) * 6 ≈ 1

G = boundary(tetrs)
Z = BEAST.ttrace(curl(X), G)

# fcr, geo = facecurrents(u, Z)
import Plotly
# Plotly.plot(patch(geo, norm.(fcr)))

Dir = Mesh(vertices(tetrs), CompScienceMeshes.celltype(G)[])
# error("stop")
Xplus = BEAST.nedelecc3d(tetrs, skeleton(tetrs,1))

bnd_tetrs = boundary(tetrs)
ttXplus = BEAST.ttrace(Xplus, bnd_tetrs)
TF, Idcs = isdivconforming(ttXplus)
@show length(Idcs)

# error("stop")
Q = BEAST.dual2forms(tetrs, skeleton(tetrs,1), Dir)



QXplus = assemble(Id, Q, Xplus)
curlX = curl(X)
QcurlX = assemble(Id, Q, curlX)

v = QXplus \ (QcurlX * u)
fcr1, geo1 = facecurrents(v, BEAST.ttrace(Xplus, bnd_tetrs));
fcr2, geo2 = facecurrents(u, BEAST.ttrace(curl(X), bnd_tetrs));

# tetrs1 = skeleton(tetrs,1)
# tetrs2 = skeleton(tetrs,2)
# ttXplus = BEAST.ttrace(Xplus, bnd_tetrs)
#
# @assert dimension(geometry(ttXplus)) == 2
# length(geometry(ttXplus)) == length(tetrs2)
#
# Conn = connectivity(tetrs1, tetrs2)
