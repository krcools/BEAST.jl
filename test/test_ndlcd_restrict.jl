using CompScienceMeshes
using BEAST
using Test
using LinearAlgebra

o, x, y, z = CompScienceMeshes.euclidianbasis(3)
tet = simplex(x,y,z,o)

rs = BEAST.NDLCDRefSpace{Float64}()
Q = BEAST.restrict(rs, tet, tet)
@test Q ≈ Matrix(1.0LinearAlgebra.I, 4, 4)
