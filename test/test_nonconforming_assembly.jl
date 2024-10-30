using BEAST, CompScienceMeshes, LinearAlgebra, Test

#= 
    This script test the strategy for the assembly of an integral operator
    for basis and test function associated with nonconforimg mesh
    It is tested by comparing the elements of interaction matrix obtained 
    using a given basis and test functions with the elements of interaction matrix obtained
    when the function are are swapped. 
=#
tol = 1e-2

m1 = meshcuboid(1.0,1.0,0.25,0.2)
m2 = meshcuboid(1.0,1.0,0.25,0.15)

X1 = BEAST.raviartthomas(m1)
X2 = BEAST.raviartthomas(m2)
t = Maxwell3D.singlelayer(wavenumber=1.0)

parent_qstrat = BEAST.DoubleNumSauterQstrat(7,8,8,8,8,8)
qstrat = BEAST.NonConformingIntegralOpQStrat(parent_qstrat)
Z21 = BEAST.assemble(t, X2, X1, quadstrat=qstrat)
Z12 = BEAST.assemble(t, X1, X2, quadstrat=qstrat)
@test maximum(norm.(Z21-transpose(Z12)))<tol

#=
# Another test scheme 

e_inc = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=1.0)
e = (n × e_inc) × n

E1 = BEAST.assemble(e, X1)
E2 = BEAST.assemble(e, X2)

#solving this overdetermined system
E_test = transpose(Z21)*E2
Z_test = transpose(Z21)*Z21
S2 = BEAST.GMRESSolver(Z1221)
j2 = solve(S2, E_test)[1]

Z11 = BEAST.assemble(t, X1, X1)
S1 = BEAST.GMRESSolver(Z11)
j1 = solve(S1, E1)[1]
=#