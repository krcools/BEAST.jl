using CompScienceMeshes
using BEAST

using Test

"""
This test is achieved by projecting the coefficients of 1st degree 
div conforming basis functions onto the 2nd degree basis functions
and then projecting back the obtained coefficients of the 2nd degree
basis functions onto the 1st degree basis functions. The resultant
operator should be an Identity operator.
"""
local m = meshrectangle(1.0,1.0,0.5)
tol = 1e-10
X2 = BEAST.raviartthomas2(m)
X1 = raviartthomas(m)
G11 = assemble(Identity(), X1, X1)
G12 = assemble(Identity(), X1, X2)
G21 = assemble(Identity(), X2, X1)
G22 = assemble(Identity(), X2, X2)
Id = I(numfunctions(X1))
@test norm(inv(Matrix(G11))*G12*inv(Matrix(G22))*G21-Id)<tol

