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

bcs = BEAST.buffachristiansen2(Faces)
rts = BEAST.raviartthomas(Faces)

G = assemble(BEAST.NCross(), bcs, rts)
