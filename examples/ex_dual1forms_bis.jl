using CompScienceMeshes
using BEAST
using LinearAlgebra

const Id = BEAST.Identity()
const Nx = BEAST.NCross()

function Base.in(mesh::CompScienceMeshes.AbstractMesh)
    cells_mesh = sort.(mesh)
    function f(cell)
        sort(cell) in cells_mesh
    end
end

function add!(fn, x, space, idcs)
    for (m,bf) in enumerate(space.fns)
        for sh in bf
            cellid = idcs[sh.cellid]
            BEAST.add!(fn, cellid, sh.refid, sh.coeff * x[m])
        end
    end
end


function compress!(space)
    T = scalartype(space)
    for (i,fn) in pairs(space.fns)
        shapes = Dict{Tuple{Int,Int},T}()
        for shape in fn
            v = get(shapes, (shape.cellid, shape.refid), zero(T))
            shapes[(shape.cellid, shape.refid)] = v + shape.coeff
        end
        space.fns[i] = [BEAST.Shape(k[1], k[2], v) for (k,v) in shapes]
    end
end

import Plotly
function showfn(space,i)
    geo = geometry(space)
    T = coordtype(geo)
    X = Dict{Int,T}()
    Y = Dict{Int,T}()
    Z = Dict{Int,T}()
    U = Dict{Int,T}()
    V = Dict{Int,T}()
    W = Dict{Int,T}()
    for sh in space.fns[i]
        chrt = chart(geo, cells(geo)[sh.cellid])
        nbd = center(chrt)
        vals = refspace(space)(nbd)
        x,y,z = cartesian(nbd)
        u,v,w = vals[sh.refid].value
        X[sh.cellid] = x
        Y[sh.cellid] = y
        Z[sh.cellid] = z
        U[sh.cellid] = get(U,sh.cellid,zero(T)) + sh.coeff * u
        V[sh.cellid] = get(V,sh.cellid,zero(T)) + sh.coeff * v
        W[sh.cellid] = get(W,sh.cellid,zero(T)) + sh.coeff * w
    end
    X = collect(values(X))
    Y = collect(values(Y))
    Z = collect(values(Z))
    U = collect(values(U))
    V = collect(values(V))
    W = collect(values(W))
    Plotly.cone(x=X,y=Y,z=Z,u=U,v=V,w=W)
end

Tetrs = CompScienceMeshes.tetmeshsphere(1.0, 0.45)
tetrs = barycentric_refinement(Tetrs)

v2t, v2n = CompScienceMeshes.vertextocellmap(tetrs)

Bnd = boundary(Tetrs)
bnd = boundary(tetrs)

Dir = Bnd
dir = bnd

Faces = submesh(!in(Bnd), CompScienceMeshes.skeleton_fast(Tetrs,2))
@show length(Faces)

F = 396
Face = cells(Faces)[F]
pos = cartesian(center(chart(Faces, Face)))

idcs1 = v2t[Face[1],1:v2n[Face[1]]]
idcs2 = v2t[Face[2],1:v2n[Face[2]]]
idcs3 = v2t[Face[3],1:v2n[Face[3]]]

supp1 = tetrs.mesh[idcs1]
supp2 = tetrs.mesh[idcs2]
supp3 = tetrs.mesh[idcs3]

supp1_edges = CompScienceMeshes.skeleton_fast(supp1, 1)
supp2_edges = CompScienceMeshes.skeleton_fast(supp2, 1)
supp3_edges = CompScienceMeshes.skeleton_fast(supp3, 1)

supp1_nodes = CompScienceMeshes.skeleton_fast(supp1, 0)
supp2_nodes = CompScienceMeshes.skeleton_fast(supp2, 0)
supp3_nodes = CompScienceMeshes.skeleton_fast(supp3, 0)

port_edges = supp1_edges
port_edges = submesh(in(supp2_edges), port_edges)
port_edges = submesh(in(supp3_edges), port_edges)
@show length(port_edges)
@assert 1 ≤ length(port_edges) ≤ 2



dir1_faces = submesh(in(dir), boundary(supp1))
dir2_faces = submesh(in(dir), boundary(supp2))
dir3_faces = submesh(in(dir), boundary(supp3))

dir1_compl = submesh(!in(dir), boundary(supp1))
dir2_compl = submesh(!in(dir), boundary(supp2))
dir3_compl = submesh(!in(dir), boundary(supp3))

dir1_compl_edges = CompScienceMeshes.skeleton_fast(dir1_compl, 1)
dir2_compl_edges = CompScienceMeshes.skeleton_fast(dir2_compl, 1)
dir3_compl_edges = CompScienceMeshes.skeleton_fast(dir3_compl, 1)

dir1_compl_nodes = CompScienceMeshes.skeleton_fast(dir1_compl, 0)
dir2_compl_nodes = CompScienceMeshes.skeleton_fast(dir2_compl, 1)
dir3_compl_nodes = CompScienceMeshes.skeleton_fast(dir3_compl, 2)

supp23 = submesh(in(boundary(supp2)), boundary(supp3))
supp31 = submesh(in(boundary(supp3)), boundary(supp1))
supp12 = submesh(in(boundary(supp1)), boundary(supp2))

supp23_edges = CompScienceMeshes.skeleton_fast(supp23, 1)
supp31_edges = CompScienceMeshes.skeleton_fast(supp31, 1)
supp12_edges = CompScienceMeshes.skeleton_fast(supp12, 1)

supp23_nodes = CompScienceMeshes.skeleton_fast(supp23, 0)
supp31_nodes = CompScienceMeshes.skeleton_fast(supp31, 0)
supp12_nodes = CompScienceMeshes.skeleton_fast(supp12, 0)

dir23_edges = submesh(in(boundary(dir2_faces)), boundary(dir3_faces))
dir31_edges = submesh(in(boundary(dir3_faces)), boundary(dir1_faces))
dir12_edges = submesh(in(boundary(dir1_faces)), boundary(dir2_faces))

dir23_compl_edges = submesh(!in(dir23_edges), boundary(supp23))
dir31_compl_edges = submesh(!in(dir31_edges), boundary(supp31))
dir12_compl_edges = submesh(!in(dir12_edges), boundary(supp12))

dir23_compl_nodes = CompScienceMeshes.skeleton_fast(dir23_compl_edges, 0)
dir31_compl_nodes = CompScienceMeshes.skeleton_fast(dir31_compl_edges, 0)
dir12_compl_nodes = CompScienceMeshes.skeleton_fast(dir12_compl_edges, 0)

supp23_int_edges = submesh(!in(dir23_compl_edges), supp23_edges)
supp31_int_edges = submesh(!in(dir31_compl_edges), supp31_edges)
supp12_int_edges = submesh(!in(dir12_compl_edges), supp12_edges)

supp23_int_nodes = submesh(!in(dir23_compl_nodes), supp23_nodes)
supp31_int_nodes = submesh(!in(dir31_compl_nodes), supp31_nodes)
supp12_int_nodes = submesh(!in(dir12_compl_nodes), supp12_nodes)

X23_prt = BEAST.nedelec(supp23, port_edges)
X31_prt = BEAST.nedelec(supp31, port_edges)
X12_prt = BEAST.nedelec(supp12, port_edges)

X23_int = BEAST.nedelec(supp23, supp23_int_edges)
X31_int = BEAST.nedelec(supp31, supp31_int_edges)
X12_int = BEAST.nedelec(supp12, supp12_int_edges)

curl_X23_prt = divergence(n × X23_prt)
curl_X31_prt = divergence(n × X31_prt)
curl_X12_prt = divergence(n × X12_prt)

curl_X23_int = divergence(n × X23_int)
curl_X31_int = divergence(n × X31_int)
curl_X12_int = divergence(n × X12_int)

Y23 = lagrangec0d1(supp23, supp23_int_nodes)
Y31 = lagrangec0d1(supp31, supp31_int_nodes)
Y12 = lagrangec0d1(supp12, supp12_int_nodes)

grad_Y23 = n × curl(Y23)
grad_Y31 = n × curl(Y31)
grad_Y12 = n × curl(Y12)

# Step 1: set port flux and extend to dual faces
x0 = ones(length(port_edges)) / length(port_edges)
for (i,edge) in enumerate(port_edges)
    tgt = tangents(center(chart(port_edges, edge)),1)
    if dot(normal(chart(Faces, Face)), tgt) < 0
        x0[i] *= -1
    end
end

A = assemble(Id, curl_X23_int, curl_X23_int)
a = -assemble(Id, curl_X23_int, curl_X23_prt) * x0
B = assemble(Id, grad_Y23, X23_int)
b = -assemble(Id, grad_Y23, X23_prt) * x0
nb = length(b)
Z = zeros(eltype(B), nb, nb)
u = [A B'; B Z] \ [a; b]
x23 = u[1:end-nb]

A = assemble(Id, curl_X31_int, curl_X31_int)
a = -assemble(Id, curl_X31_int, curl_X31_prt) * x0
B = assemble(Id, grad_Y31, X31_int)
b = -assemble(Id, grad_Y31, X31_prt) * x0
nb = length(b)
Z = zeros(eltype(B), nb, nb)
u = [A B'; B Z] \ [a; b]
x31 = u[1:end-nb]

A = assemble(Id, curl_X12_int, curl_X12_int)
a = -assemble(Id, curl_X12_int, curl_X12_prt) * x0
B = assemble(Id, grad_Y12, X12_int)
b = -assemble(Id, grad_Y12, X12_prt) * x0
nb = length(b)
Z = zeros(eltype(B), nb, nb)
u = [A B'; B Z] \ [a; b]
x12 = u[1:end-nb]


# Step 2: extend to the volume
supp1_int_edges = submesh(!in(dir1_compl_edges), supp1_edges)
supp2_int_edges = submesh(!in(dir2_compl_edges), supp2_edges)
supp3_int_edges = submesh(!in(dir3_compl_edges), supp3_edges)

supp1_int_nodes = submesh(!in(dir1_compl_nodes), supp1_nodes)
supp2_int_nodes = submesh(!in(dir2_compl_nodes), supp2_nodes)
supp3_int_nodes = submesh(!in(dir3_compl_nodes), supp3_nodes)

supp1_drv_edges = CompScienceMeshes.union(port_edges, supp31_int_edges)
supp1_drv_edges = CompScienceMeshes.union(supp1_drv_edges, supp12_int_edges)
supp2_drv_edges = CompScienceMeshes.union(port_edges, supp12_int_edges)
supp2_drv_edges = CompScienceMeshes.union(supp2_drv_edges, supp23_int_edges)
supp3_drv_edges = CompScienceMeshes.union(port_edges, supp23_int_edges)
supp3_drv_edges = CompScienceMeshes.union(supp3_drv_edges, supp31_int_edges)

X1_drv = BEAST.nedelecc3d(supp1, supp1_drv_edges)
X2_drv = BEAST.nedelecc3d(supp2, supp2_drv_edges)
X3_drv = BEAST.nedelecc3d(supp3, supp3_drv_edges)

X1_int = BEAST.nedelecc3d(supp1, supp1_int_edges)
X2_int = BEAST.nedelecc3d(supp2, supp2_int_edges)
X3_int = BEAST.nedelecc3d(supp3, supp3_int_edges)

curl_X1_drv = curl(X1_drv)
curl_X2_drv = curl(X2_drv)
curl_X3_drv = curl(X3_drv)

curl_X1_int = curl(X1_int)
curl_X2_int = curl(X2_int)
curl_X3_int = curl(X3_int)

Y1 = BEAST.lagrangec0d1(supp1, supp1_int_nodes)
Y2 = BEAST.lagrangec0d1(supp2, supp2_int_nodes)
Y3 = BEAST.lagrangec0d1(supp3, supp3_int_nodes)

grad_Y1 = BEAST.gradient(Y1)
grad_Y2 = BEAST.gradient(Y2)
grad_Y3 = BEAST.gradient(Y3)

A = assemble(Id, curl_X1_int, curl_X1_int)
a = -assemble(Id, curl_X1_int, curl_X1_drv) * [x0; x31; x12]
B = assemble(Id, grad_Y1, X1_int)
b = -assemble(Id, grad_Y1, X1_drv) * [x0; x31; x12]
nb = length(b)
Z = zeros(eltype(B), nb, nb)
u = [A B'; B Z] \ [a; b]
x1 = u[1:end-nb]

A = assemble(Id, curl_X2_int, curl_X2_int)
a = -assemble(Id, curl_X2_int, curl_X2_drv) * [x0; x12; x23]
B = assemble(Id, grad_Y2, X2_int)
b = -assemble(Id, grad_Y2, X2_drv) * [x0; x12; x23]
nb = length(b)
Z = zeros(eltype(B), nb, nb)
u = [A B'; B Z] \ [a; b]
x2 = u[1:end-nb]

A = assemble(Id, curl_X3_int, curl_X3_int)
a = -assemble(Id, curl_X3_int, curl_X3_drv) * [x0; x23; x31]
B = assemble(Id, grad_Y3, X3_int)
b = -assemble(Id, grad_Y3, X3_drv) * [x0; x23; x31]
nb = length(b)
Z = zeros(eltype(B), nb, nb)
u = [A B'; B Z] \ [a; b]
x3 = u[1:end-nb]


# inject in the global space
fn = BEAST.Shape{Float64}[]
add!(fn, [x0; x31; x12], X1_drv, idcs1)
add!(fn, x1, X1_int, idcs1)

add!(fn, [x0; x12; x23], X2_drv, idcs2)
add!(fn, x2, X2_int, idcs2)

add!(fn, [x0; x23; x31], X3_drv, idcs3)
add!(fn, x3, X3_int, idcs3)

pos = cartesian(CompScienceMeshes.center(chart(Faces, Face)))
space = BEAST.NDLCCBasis(tetrs, [fn], [pos])
