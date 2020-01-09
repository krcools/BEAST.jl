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

Tetrs = CompScienceMeshes.tetmeshsphere(1.0, 0.35)
tetrs = barycentric_refinement(Tetrs)

cells_Bndry = [sort(c) for c in cells(skeleton(boundary(Tetrs),1))]
Edges = submesh(skeleton(Tetrs,1)) do Edge
    sort(Edge) in cells_Bndry ? false : true
end
Edges = Mesh(vertices(Edges), cells(Edges))
@show numcells(Edges)

# pred = CompScienceMeshes.interior_tpredicate(Tetrs)
# AllFaces = skeleton(Tetrs,2)
# sm = submesh(pred, AllFaces)

E = 102

Edge = cells(Edges)[E]
pos = cartesian(center(chart(Edges, Edge)))
# v = argmin(norm.(vertices(tetrs) .- Ref(pos)))
support1 = submesh(tetr -> Edge[1] in tetr, tetrs.mesh)
support2 = submesh(tetr -> Edge[2] in tetr, tetrs.mesh)
# port = submesh(face -> v in face, boundary(support1))
# port2 = submesh(face -> v in face, boundary(support2))
port = submesh(face -> sort(face) in sort.(boundary(support2)), boundary(support1))
for p in port
    @assert sort(p) in sort.(port2)
end
@show length(port)

support = CompScienceMeshes.union(support1, support2)
@assert CompScienceMeshes.isoriented(support)
@show length(support)
@assert length(skeleton(support1,2)) + length(skeleton(support2,2)) ==
    length(skeleton(support,2)) + length(port)
edges = skeleton(support,1)

cells_bnd_edges = [sort(c) for c in skeleton(boundary(support),1)]
cells_prt_edges = [sort(c) for c in skeleton(port,1)]
cells_dir_edges = [sort(c) for c in skeleton(boundary(tetrs),1)]

bnd_support = boundary(support)
cells_bnd_tetrs = sort.(boundary(tetrs))
dir_cap_bnd = submesh(face -> sort(face) in cells_bnd_tetrs, bnd_support)
cells_bnd_dir_cap_bnd = sort.(boundary(dir_cap_bnd))

int_edges = submesh(edges) do edge
    sort(edge) in cells_bnd_dir_cap_bnd && return false
    sort(edge) in cells_dir_edges && return true
    sort(edge) in cells_bnd_edges && return false
    sort(edge) in cells_prt_edges && return false
    return true
end
@show length(int_edges)

Nd_int = BEAST.nedelecc3d(support, int_edges)
Q = assemble(Id, Nd_int, Nd_int)

Q2, store = BEAST.allocatestorage(Id, Nd_int, Nd_int,
    Val{:bandedstorage}, BEAST.LongDelays{:ignore})
BEAST.assemble_local_mixed!(Id, Nd_int, Nd_int, store)

@assert Q ≈ Q2

Q3 = assemble(Id, curl(Nd_int), curl(Nd_int))
rank(Q3)

cells_bnd_faces = [sort(c) for c in boundary(support)]
cells_prt_faces = [sort(c) for c in port]
int_faces = submesh(skeleton(support,2)) do face
    sort(face) in cells_bnd_tetrs && return true
    sort(face) in cells_bnd_faces && return false
    sort(face) in cells_prt_faces && return false
    return true
end
@show length(int_faces)
@assert length(int_faces) + length(boundary(support)) + length(port) ==
    length(skeleton(support,2)) + length(dir_cap_bnd)

RT_int = BEAST.nedelecd3d(support, int_faces)
for (m,fn) in enumerate(RT_int.fns)
    ct = count(sh -> sh.cellid <= length(support1), fn)
    # @show m, ct
    @assert ct == 2 || ct == 0
end
Q4 = assemble(Id, divergence(RT_int), divergence(RT_int))
numfunctions(RT_int) - rank(Q4)

RT_prt = BEAST.nedelecd3d(support, port)
for fn in RT_prt.fns
    @assert count(sh -> sh.cellid <= length(support1), fn) == 1
    @assert count(sh -> sh.cellid > length(support1), fn) == 1
end
Q5 = assemble(Id, divergence(RT_int), divergence(RT_int))
x0 = ones(length(port)) / length(port)
tgt = vertices(Edges)[Edge[1]] - vertices(Edges)[Edge[2]]
for (i,face) in enumerate(port)
    x0[i] *= sign(dot(normal(chart(port,face)), tgt))
    # if RT_prt.fns[i][1].coeff < 0
    #     x0[i] *= -1
    # end
end
Q6 = assemble(Id, divergence(RT_int), divergence(RT_prt))
d = -Q6 * x0
x1 = pinv(Q5) * d

# L0 = lagrangecxd0(support)
# Q7 = assemble(Id, L0, divergence(RT_int))
# Q8 = assemble(Id, L0, L0)
# ch1 = []
# ch2 = []
# ch = [ch1; ch2]
# x1 = pinv(Q7) * ()

fn = BEAST.Shape{Float64}[]
add!(fn, x0, RT_prt)
add!(fn, x1, RT_int)
Y1 = BEAST.NDLCDBasis(support, [fn], [pos])
divY1 = divergence(Y1); compress!(divY1)


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

error("stop")
Y = BEAST.dual2forms(Tetrs, Edges)
