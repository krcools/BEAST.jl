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

function extend_edge_to_face(supp, dirichlet, x_prt, port_edges)
    bnd_supp = boundary(supp)
    supp_edges = CompScienceMeshes.skeleton_fast(supp, 1)
    supp_nodes = CompScienceMeshes.skeleton_fast(supp, 0)

    dir_compl_edges = submesh(!in(dirichlet), bnd_supp)
    dir_compl_nodes = CompScienceMeshes.skeleton_fast(dir_compl_edges, 0)

    int_edges = submesh(!in(dir_compl_edges), supp_edges)
    int_nodes = submesh(!in(dir_compl_nodes), supp_nodes)

    Nd_prt = BEAST.nedelec(supp, port_edges)
    Nd_int = BEAST.nedelec(supp, int_edges)
    Lg_int = BEAST.lagrangec0d1(supp, int_nodes)

    curl_Nd_prt = divergence(n × Nd_prt)
    curl_Nd_int = divergence(n × Nd_int)
    grad_Lg_int = n × curl(Lg_int)

    A = assemble(Id, curl_Nd_int, curl_Nd_int)
    B = assemble(Id, grad_Lg_int, Nd_int)

    a = -assemble(Id, curl_Nd_int, curl_Nd_prt) * x_prt
    b = -assemble(Id, grad_Lg_int, Nd_prt) * x_prt

    Z = zeros(eltype(b), length(b), length(b))
    u = [A B'; B Z] \ [a;b]
    x_int = u[1:end-length(b)]

    return x_int, int_edges, Nd_int
end



function extend_face_to_tetr(supp, dirichlet, x_prt, port_edges)
    bnd_supp = boundary(supp)
    supp_edges = CompScienceMeshes.skeleton_fast(supp, 1)
    supp_nodes = CompScienceMeshes.skeleton_fast(supp, 0)

    dir_compl_faces = submesh(!in(dirichlet), bnd_supp)
    dir_compl_edges = CompScienceMeshes.skeleton_fast(dir_compl_faces, 1)
    dir_compl_nodes = CompScienceMeshes.skeleton_fast(dir_compl_faces, 0)

    int_edges = submesh(!in(dir_compl_edges), supp_edges)
    int_nodes = submesh(!in(dir_compl_nodes), supp_nodes)

    Nd_prt = BEAST.nedelecc3d(supp, port_edges)
    Nd_int = BEAST.nedelecc3d(supp, int_edges)
    Lg_int = BEAST.lagrangec0d1(supp, int_nodes)

    curl_Nd_prt = curl(Nd_prt)
    curl_Nd_int = curl(Nd_int)
    grad_Lg_int = BEAST.gradient(Lg_int)

    A = assemble(Id, curl_Nd_int, curl_Nd_int)
    B = assemble(Id, grad_Lg_int, Nd_int)

    a = -assemble(Id, curl_Nd_int, curl_Nd_prt) * x_prt
    b = -assemble(Id, grad_Lg_int, Nd_prt) * x_prt

    Z = zeros(eltype(b), length(b), length(b))
    u = [A B'; B Z] \ [a;b]
    x_int = u[1:end-length(b)]

    return x_int, int_edges, Nd_int
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

dir1_faces = submesh(in(dir), boundary(supp1))
dir2_faces = submesh(in(dir), boundary(supp2))
dir3_faces = submesh(in(dir), boundary(supp3))

supp23 = submesh(in(boundary(supp2)), boundary(supp3))
supp31 = submesh(in(boundary(supp3)), boundary(supp1))
supp12 = submesh(in(boundary(supp1)), boundary(supp2))

dir23_edges = submesh(in(boundary(dir2_faces)), boundary(dir3_faces))
dir31_edges = submesh(in(boundary(dir3_faces)), boundary(dir1_faces))
dir12_edges = submesh(in(boundary(dir1_faces)), boundary(dir2_faces))
port_edges = boundary(supp23)
port_edges = submesh(in(boundary(supp31)), port_edges)
port_edges = submesh(in(boundary(supp12)), port_edges)
@assert 1 ≤ length(port_edges) ≤ 2

# Step 1: set port flux and extend to dual faces
x0 = ones(length(port_edges)) / length(port_edges)
for (i,edge) in enumerate(port_edges)
    tgt = tangents(center(chart(port_edges, edge)),1)
    if dot(normal(chart(Faces, Face)), tgt) < 0
        x0[i] *= -1
    end
end

x23, supp23_int_edges = extend_edge_to_face(supp23, dir23_edges, x0, port_edges)
x31, supp31_int_edges = extend_edge_to_face(supp31, dir31_edges, x0, port_edges)
x12, supp12_int_edges = extend_edge_to_face(supp12, dir12_edges, x0, port_edges)

port1_edges = CompScienceMeshes.union(port_edges, supp31_int_edges)
port1_edges = CompScienceMeshes.union(port1_edges, supp12_int_edges)

port2_edges = CompScienceMeshes.union(port_edges, supp12_int_edges)
port2_edges = CompScienceMeshes.union(port2_edges, supp23_int_edges)

port3_edges = CompScienceMeshes.union(port_edges, supp23_int_edges)
port3_edges = CompScienceMeshes.union(port3_edges, supp31_int_edges)

x1_prt = [x0; x31; x12]
x2_prt = [x0; x12; x23]
x3_prt = [x0; x23; x31]

Nd1_prt = BEAST.nedelecc3d(supp1, port1_edges)
Nd2_prt = BEAST.nedelecc3d(supp2, port2_edges)
Nd3_prt = BEAST.nedelecc3d(supp3, port3_edges)

x1_int, _, Nd1_int = extend_face_to_tetr(supp1, dir1_faces, x1_prt, port1_edges)
x2_int, _, Nd2_int = extend_face_to_tetr(supp2, dir2_faces, x2_prt, port2_edges)
x3_int, _, Nd3_int = extend_face_to_tetr(supp3, dir3_faces, x3_prt, port3_edges)

# inject in the global space
fn = BEAST.Shape{Float64}[]
add!(fn, x1_prt, Nd1_prt, idcs1)
add!(fn, x1_int, Nd1_int, idcs1)

add!(fn, x2_prt, Nd2_prt, idcs2)
add!(fn, x2_int, Nd2_int, idcs2)

add!(fn, x3_prt, Nd3_prt, idcs3)
add!(fn, x3_int, Nd3_int, idcs3)

pos = cartesian(CompScienceMeshes.center(chart(Faces, Face)))
space = BEAST.NDLCCBasis(tetrs, [fn], [pos])


tetrs, bnd, dir, v2t, v2n = BEAST.dual1forms_init(Tetrs, Dir)
D1 = BEAST.dual1forms_body(Faces, tetrs, bnd, dir, v2t, v2n)
P2 = BEAST.nedelecd3d(Tetrs, Faces)
# S = BEAST.dual1forms(Tetrs, Faces[collect(396:405)], Dir)

Edges = CompScienceMeshes.skeleton_fast(Tetrs, 1)
int_Edges = submesh(!in(skeleton(boundary(Tetrs),1)), Edges)
P1 = BEAST.nedelecc3d(Tetrs, int_Edges)

curl_D1 = curl(D1)
curl_P1 = curl(P1)
div_P2 = divergence(P2)

G = assemble(Id, curl_D1, curl_D1)
GP = assemble(Id, div_P2, div_P2)
M = assemble(Id, D1, P2)

CDP = assemble(Id, D1, curl_P1)
CPD = CDP'
GDD = assemble(Id, D1, D1)

iM = inv(M)
A1 = CPD * inv(GDD) * CDP;
A2 = assemble(Id, curl_P1, curl_P1)
A3 = CPD * iM' * inv(GDD) * iM * CDP;
