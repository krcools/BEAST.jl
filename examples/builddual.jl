using CompScienceMeshes
using BEAST

Faces = meshsphere(1.0, 0.35)
Edges = skeleton(Faces,1)

faces = barycentric_refinement(Faces)
edges = skeleton(faces,1)
# verts = skeleton(faces,0)

E = 1
Edge = cells(Edges)[1]

port_idx = numvertices(Faces) + E
ptch_idx = Edge[1]

# patch = Mesh(vertices(faces), filter(c -> ptch_idx in c, cells(faces)))
patch_idcs = Int[]
for (i,face) in enumerate(cells(faces))
    if ptch_idx in face
        push!(patch_idcs, i)
    end
end
patch = Mesh(vertices(faces), cells(faces)[patch_idcs])

port = Mesh(vertices(edges), filter(c -> port_idx in c, cells(boundary(patch))))

@show numcells(patch)
@show numcells(port)

# D, C, d, c, d0, d1,
RT_int, RT_prt, x_int, x_prt = BEAST.buildhalfbc2(patch, port, nothing)

BF = BEAST.Shape{Float64}[]
for (m,bf) in enumerate(RT_int.fns)
    for sh in bf
        cellid = patch_idcs[sh.cellid]
        BEAST.add!(BF,cellid, sh.refid, sh.coeff * x_int[m])
    end
end

for (m,bf) in enumerate(RT_prt.fns)
    for sh in bf
        cellid = patch_idcs[sh.cellid]
        BEAST.add!(BF,cellid, sh.refid, sh.coeff * x_prt[m])
    end
end


# RT_prt = raviartthomas(patch, cellpairs(patch, port))
# Lx = lagrangecxd0(patch)
# @assert numfunctions(Lx) == numcells(patch)
#
# Id = BEAST.Identity()
# div_RT_prt = divergence(RT_prt)
# Z, store = BEAST.allocatestorage(Id, Lx, div_RT_prt, Val{:bandedstorage}, BEAST.LongDelays{:ignore})
# BEAST.assemble_local_mixed!(Id, Lx, div_RT_prt, store)

bcs3 = BEAST.buffachristiansen3(Faces)
bcs2 = BEAST.buffachristiansen2(Faces)
bcs1 = BEAST.buffachristiansen(Faces)
rts = BEAST.raviartthomas(Faces)

G1 = assemble(BEAST.NCross(), bcs1, rts)
Q1 = assemble(BEAST.Identity(), divergence(bcs1), divergence(bcs1))

G2 = assemble(BEAST.NCross(), bcs2, rts)
Q2 = assemble(BEAST.Identity(), divergence(bcs2), divergence(bcs2))

G3 = assemble(BEAST.NCross(), bcs3, rts)
Q3 = assemble(BEAST.Identity(), divergence(bcs3), divergence(bcs3))

using LinearAlgebra
@show cond(G3)

using LinearAlgebra
@assert (numcells(Faces) - 1) == (size(Q,1) - rank(Q))
@assert cond(G) < 3.5

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

getch(geo,i) = chart(geo, cells(geo)[i])


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
        @show vals[sh.refid].value
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
