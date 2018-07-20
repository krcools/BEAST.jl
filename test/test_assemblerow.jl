using CompScienceMeshes, BEAST
using Test

fn = joinpath(dirname(@__FILE__),"assets","sphere35.in")
m = readmesh(fn)
t = Maxwell3D.singlelayer(wavenumber=1.0)
X = raviartthomas(m)
numfunctions(X)

##
X1 = subset(X,1:1)
numfunctions(X1)

T1 = assemble(t,X1,X)
T2 = BEAST.assemblerow(t,X1,X)
#
# T3 = assemble(t,X,X1)
# T4 = BEAST.assemblecol(t,X,X1)
#
# @test T1 == T2

# T2 = BEAST.assembleblock(t,X,X)
@test T1 == T2


I = [3,2,7]
X1 = subset(X,I)
T2 = zeros(scalartype(t,X1,X1),numfunctions(X1),numfunctions(X1))
store(v,m,n) = (T2[m,n] += v)

test_elements, test_assembly_data,
	trial_elements, trial_assembly_data,
	quadrature_data, zlocal = BEAST.assembleblock_primer(t,X,X)

BEAST.assembleblock_body!(t,
	X, I, test_elements,  test_assembly_data,
    X, I, trial_elements, trial_assembly_data,
    quadrature_data, zlocal, store)

T1 = assemble(t,X1,X1)
@test T1 == T2

# @time BEAST.assembleblock_body!(t,
# 	X, I, test_elements,  test_assembly_data,
#     X, I, trial_elements, trial_assembly_data,
#     quadrature_data, zlocal, store)
# @time assemble(t,X1,X1)

T3 = zeros(scalartype(t,X1,X1),numfunctions(X1),numfunctions(X1))
store3(v,m,n) = (T3[m,n] += v)

blkasm = BEAST.blockassembler(t,X,X)
blkasm(I,I,store3)

T4 = assemble(t,X,X)
@test T3 == T4[I,I]

# @time blkasm(I,I)
# @time assemble(t,X1,X1)
