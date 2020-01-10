using CompScienceMeshes
using BEAST
using LinearAlgebra

const Id = BEAST.Identity()

function add!(fn, x, space)
    for (m,bf) in enumerate(space.fns)
        for sh in bf
            BEAST.add!(fn, sh.cellid, sh.refid, sh.coeff * x[m])
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
            # set!(shapes, (shape.cellid, shape.refid), v + shape.coeff)
        end
        space.fns[i] = [BEAST.Shape(k[1], k[2], v) for (k,v) in shapes]
    end
end

import PlotlyJS
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
        # @show vals[sh.refid].value
        u,v,w = vals[sh.refid].value
        # @show x, y, z
        # @show u, v, w
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
    PlotlyJS.cone(x=X,y=Y,z=Z,u=U,v=V,w=W)
end

Tetrs = CompScienceMeshes.tetmeshsphere(1.0, 0.35)
tetrs = barycentric_refinement(Tetrs)

bnd_Tetrs = boundary(Tetrs)
bnd_tetrs = boundary(tetrs)
srt_bnd_Tetrs = sort.(bnd_Tetrs)
srt_bnd_tetrs = sort.(bnd_tetrs)

Faces = submesh(Face -> !(sort(Face) in srt_bnd_Tetrs), skeleton(Tetrs,2))
@show length(Faces)

F = 603
Face = cells(Faces)[F]
pos = cartesian(center(chart(Faces, Face)))

supp1 = submesh(tet -> Face[1] in tet, tetrs.mesh)
supp2 = submesh(tet -> Face[2] in tet, tetrs.mesh)
supp3 = submesh(tet -> Face[3] in tet, tetrs.mesh)

supp = CompScienceMeshes.union(supp1, supp2)
supp = CompScienceMeshes.union(supp, supp3)
# PlotlyJS.plot(patch(boundary(supp)))

port = skeleton(supp1, 1)
port = submesh(edge -> sort(edge) in sort.(skeleton(supp2,1)), port)
port = submesh(edge -> sort(edge) in sort.(skeleton(supp3,1)), port)
@show length(port)
@assert 1 ≤ length(port) ≤ 2

dir = submesh(face -> sort(face) in srt_bnd_tetrs, boundary(supp))
rim = boundary(dir)
@show length(dir)

srt_port = sort.(port)
srt_dir_edges = sort.(skeleton(dir,1))
srt_bnd_edges = sort.(skeleton(boundary(supp),1))
srt_rim = sort.(rim)
int_edges = submesh(skeleton(supp,1)) do edge
    sort(edge) in srt_port && return false
    sort(edge) in srt_rim && return false
    sort(edge) in srt_dir_edges && return true
    sort(edge) in srt_bnd_edges && return false
    return true
end
@show length(int_edges)
@assert length(skeleton(supp,0)) - length(skeleton(supp,1)) +
    length(skeleton(supp,2)) - length(supp) == 1

Nd_prt = BEAST.nedelecc3d(supp, port)
Nd_int = BEAST.nedelecc3d(supp, int_edges)
@show numfunctions(Nd_int)

x0 = ones(length(port)) / length(port)
for (i,edge) in enumerate(port)
    tgt = tangents(center(chart(port,edge)),1)
    if dot(normal(chart(Faces,Face)), tgt) < 0
        x0[i] *= -1
    end
end

A = assemble(Id, curl(Nd_int), curl(Nd_int))
N = nullspace(A)
@show size(N,2)
a = -assemble(Id, curl(Nd_int), curl(Nd_prt)) * x0
x1 = pinv(A) * a
@show norm(N'*a)
@show norm(A*x1-a)

srt_bnd_verts = sort.(skeleton(boundary(supp),0))
srt_prt_verts = sort.(skeleton(port,0))
int_verts = submesh(skeleton(supp,0)) do vert
    sort(vert) in srt_bnd_verts && return false
    sort(vert) in srt_prt_verts && return false
    return true
end
@show length(int_verts)

L0_int = lagrangec0d1(supp, int_verts)
@show numfunctions(L0_int)
grad_L0_int = BEAST.gradient(L0_int)
@assert numfunctions(grad_L0_int) == numfunctions(L0_int)
B = assemble(Id, grad_L0_int, Nd_int)

B2, store = BEAST.allocatestorage(Id, grad_L0_int, Nd_int,
    Val{:bandedstorage}, BEAST.LongDelays{:ignore})
BEAST.assemble_local_mixed!(Id, grad_L0_int, Nd_int, store)
@assert isapprox(B, B2, atol=1e-8)
@show rank(B2)

b = -assemble(Id, grad_L0_int, Nd_prt) * x0
p = (B*N) \ (b-B*x1)
x1 = x1 + N*p

@show norm(A*x1-a)
@show norm(B*x1-b)

fn = BEAST.Shape{Float64}[]
add!(fn, x0, Nd_prt)
add!(fn, x1, Nd_int)

Y1 = BEAST.NDLCCBasis(supp, [fn],[pos])
curlY1 = curl(Y1); compress!(curl(Y1));

include(joinpath(@__DIR__, "utils/edge_values.jl"))
EV = edge_values(Y1,1)
check_edge_values(EV)

Y = BEAST.dual1forms(Tetrs, Faces)
