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

Tetrs = CompScienceMeshes.tetmeshsphere(1.0, 0.15)
tetrs = barycentric_refinement(Tetrs)

cells_Bndry = [sort(c) for c in cells(skeleton(boundary(Tetrs),1))]

Bnd = boundary(Tetrs)
Neu = Bnd
Dir = Mesh(vertices(Bnd), CompScienceMeshes.celltype(Bnd)[])

bnd = boundary(tetrs)
neu = bnd
dir = Mesh(vertices(bnd), CompScienceMeshes.celltype(bnd)[])

# srt_Dir = sort.(Dir)
# Edges = submesh(skeleton(Tetrs,1)) do Edge
#     sort(Edge) in srt_Dir && return false
#     return true
# end

Edges = submesh(!in(Dir), skeleton(Tetrs,1))

# Edges = skeleton(Tetrs,1)
# srt_bnd_Faces = sort.(boundary(Tetrs))
# srt_neu = sort.(neu)
# Faces = submesh(skeleton(Tetrs,2)) do Face
#     !(sort(Face) in srt_bnd_Faces)
# end

# srt_bnd_Nodes = sort.(skeleton(boundary(Tetrs),0))
# Nodes = submesh(skeleton(Tetrs,0)) do node
#     !(sort(node) in srt_bnd_Nodes)
# end

Nodes = submesh(!in(skeleton(boundary(Tetrs),0)), skeleton(Tetrs,0))

@show numcells(Edges)
# @show length(Faces)

# pred = CompScienceMeshes.interior_tpredicate(Tetrs)
# AllFaces = skeleton(Tetrs,2)
# sm = submesh(pred, AllFaces)

E = 1

Edge = cells(Edges)[E]
pos = cartesian(CompScienceMeshes.center(chart(Edges, E)))
# v = argmin(norm.(vertices(tetrs) .- Ref(pos)))
support1 = submesh((m,tetr) -> Edge[1] in CompScienceMeshes.indices(m,tetr), tetrs.mesh)
support2 = submesh((m,tetr) -> Edge[2] in CompScienceMeshes.indices(m,tetr), tetrs.mesh)
# port = submesh(face -> v in face, boundary(support1))
# port2 = submesh(face -> v in face, boundary(support2))
# port = submesh(face -> sort(face) in sort.(boundary(support2)), boundary(support1))
port = submesh(in(boundary(support2)), boundary(support1))
# for p in port
#     @assert sort(p) in sort.(port2)
# end
@show length(port)

support = CompScienceMeshes.union(support1, support2)
@assert CompScienceMeshes.isoriented(support)
@show length(support)
@assert length(skeleton(support1,2)) + length(skeleton(support2,2)) ==
    length(skeleton(support,2)) + length(port)
edges = skeleton(support,1)

# cells_bnd_edges = [sort(c) for c in skeleton(boundary(support),1)]
# cells_prt_edges = [sort(c) for c in skeleton(port,1)]
# cells_dir_edges = [sort(c) for c in skeleton(boundary(tetrs),1)]

bnd_support = boundary(support)
# cells_bnd_tetrs = sort.(boundary(tetrs))

# srt_dir = sort.(dir)
# dir_cap_bnd = submesh(face -> sort(face) in sort.(dir), bnd_support)
dir_cap_bnd = submesh(in(dir), bnd_support)
# cells_bnd_dir_cap_bnd = sort.(boundary(dir_cap_bnd))

in_bnd_dir_cap_bnd = in(boundary(dir_cap_bnd))
in_dir_edges = in(skeleton(boundary(tetrs),1))
in_bnd_edges = in(skeleton(boundary(support),1))
in_prt_edges = in(skeleton(port,1))

int_edges = submesh(edges) do m,edge
    in_bnd_dir_cap_bnd(m,edge) && return false
    in_bnd_edges(m,edge) && return true
    in_dir_edges(m,edge) && return false
    in_prt_edges(m,edge) && return false
    return true
end

# int_edges = submesh(edges) do edge
#     sort(edge) in cells_bnd_dir_cap_bnd && return false
#     sort(edge) in cells_dir_edges && return true
#     sort(edge) in cells_bnd_edges && return false
#     sort(edge) in cells_prt_edges && return false
#     return true
# end
@show length(int_edges)

Nd_int = BEAST.nedelecc3d(support, int_edges)
Q = assemble(Id, Nd_int, Nd_int)

freeze, store = BEAST.allocatestorage(Id, Nd_int, Nd_int,
    Val{:bandedstorage}, BEAST.LongDelays{:ignore})
BEAST.assemble_local_mixed!(Id, Nd_int, Nd_int, store)
Q2 = freeze()

@assert Q ≈ Q2

Q3 = assemble(Id, curl(Nd_int), curl(Nd_int))
rank(Q3)

# cells_bnd_faces = [sort(c) for c in boundary(support)]
# cells_prt_faces = [sort(c) for c in port]
# srt_dir_cap_bnd = sort.(dir_cap_bnd)

in_dir_cap_bnd = in(dir_cap_bnd)
in_bnd_faces = in(boundary(support))
in_prt_faces = in(port)

int_faces = submesh(skeleton(support,2)) do m,face
    in_dir_cap_bnd(m,face) && return true
    in_bnd_faces(m,face) && return false
    in_prt_faces(m,face) && return false
    return true
end

# int_faces = submesh(skeleton(support,2)) do face
#     # sort(face) in cells_bnd_tetrs && return true
#     sort(face) in srt_dir_cap_bnd && return true
#     sort(face) in cells_bnd_faces && return false
#     sort(face) in cells_prt_faces && return false
#     return true
# end
@show length(int_faces)
@assert length(int_faces) + length(boundary(support)) + length(port) ==
    length(skeleton(support,2)) + length(dir_cap_bnd)

RT_int = BEAST.nedelecd3d(support, int_faces)
for (m,fn) in enumerate(RT_int.fns)
    ct = count(sh -> sh.cellid <= length(support1), fn)
    # @show m, ct
    @assert 0 ≤ ct ≤ 2
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
x1 = pinv(Matrix(Q5)) * d

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
        nbd = CompScienceMeshes.center(chrt)
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
    Plotly.cone(x=X,y=Y,z=Z,u=U,v=V,w=W)
end

# error("stop")
Dir =  Mesh(vertices(Tetrs), CompScienceMeshes.celltype(int_faces)[])
# error()

tetrs, bnd, dir, v2t, v2n = BEAST.dualforms_init(Tetrs, Dir)
Y = BEAST.dual2forms_body(Edges[collect(1:10)],  tetrs, bnd, dir, v2t, v2n)

# Y = BEAST.dual2forms(Tetrs, Edges, Dir)
nothing
