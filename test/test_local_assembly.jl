@info "Executing test_local_assembly.jl"
using CompScienceMeshes
using BEAST

using Test

m = meshrectangle(1.0, 1.0, 0.5, 3)
nodes = skeleton(m,0)
int_nodes = submesh(!in(skeleton(boundary(m),0)), nodes)
# srt_bnd_nodes = sort.(skeleton(boundary(m),0))
# int_nodes = submesh(nodes) do node
#     !(sort(node) in srt_bnd_nodes)
# end

@test length(int_nodes) == 1

L0 = BEAST.lagrangec0d1(m, int_nodes)
@test numfunctions(L0) == 1
@test length(L0.fns[1]) == 6
@test length(m) == 8

Lx = BEAST.lagrangecxd0(m)
@test numfunctions(Lx) == 8

Id = BEAST.Identity()
fr1, st1 = BEAST.allocatestorage(Id, L0, Lx, Val{:bandedstorage}, BEAST.LongDelays{:ignore})
BEAST.assemble_local_matched!(Id, L0, Lx, st1)
Q1 = fr1()

fr2, st2 = BEAST.allocatestorage(Id, L0, Lx, Val{:bandedstorage}, BEAST.LongDelays{:ignore})
BEAST.assemble_local_mixed!(Id, L0, Lx, st2)
Q2 = fr2()

@test isapprox(Q1, Q2, atol=1e-8)


RT = raviartthomas(m)
BC = buffachristiansen(m)

fr1, st1 = BEAST.allocatestorage(Id, BC, RT, Val{:bandedstorage}, BEAST.LongDelays{:ignore})
BEAST.assemble_local_refines!(Id, BC, RT, st1)
Q1 = fr1()

fr2, st2 = BEAST.allocatestorage(Id, BC, RT, Val{:bandedstorage}, BEAST.LongDelays{:ignore})
BEAST.assemble_local_mixed!(Id, BC, RT, st2)
Q2 = fr2()
@test isapprox(Q1, Q2, atol=1e-8)
