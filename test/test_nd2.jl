using CompScienceMeshes
using BEAST
using LinearAlgebra

using Test


#= """
This test is achieved by projecting the coefficients of 1st degree 
curl conforming basis functions onto the 2nd degree basis functions
and then projecting back the obtained coefficients of the 2nd degree
basis functions onto the 1st degree basis functions. The resultant
operator should be an Identity operator.
""" =#

m = meshrectangle(1.0,1.0,0.5)
tol = 1e-10
X2 = BEAST.nedelec2(m)
X1 = BEAST.nedelec(m)
G11 = assemble(Identity(), X1, X1)
G12 = assemble(Identity(), X1, X2)
G21 = assemble(Identity(), X2, X1)
G22 = assemble(Identity(), X2, X2)
Id = Matrix{eltype(G11)}(LinearAlgebra.I, numfunctions(X1), numfunctions(X1))
@test norm(inv(Matrix(G11))*G12*inv(Matrix(G22))*G21-Id)<tol
