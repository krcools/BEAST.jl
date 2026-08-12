using Test
using CompScienceMeshes
import BEAST
import BEAST.BlockArrays

fn = joinpath(dirname(pathof(BEAST)), "../examples/assets/sphere45.in")

m1 = readmesh(fn)
m2 = m1[Int[]]

X = BEAST.DirectProductSpace([
    raviartthomas(m) for m in [m1, m2]
])

T = Maxwell3D.singlelayer(gamma=1.0)

@hilbertspace j[1:2]
@hilbertspace k[1:2]

a = (
    T[k[1],j[1]] +
    T[k[1],j[2]] +
    T[k[2],j[2]] +
    T[k[2],j[1]]
)

A = assemble(a, X, X)

M = AbstractMatrix(A)

n1 = numfunctions(X[1])
n2 = numfunctions(X[2])

@test n2 == 0

@test BlockArrays.blocksize(M) == (2,2)

@test BlockArrays.blocksizes(M) ==
    [(n1,n1) (n1,n2);
        (n2,n1) (n2,n2)]

