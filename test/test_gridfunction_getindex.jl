using Test
using BEAST
using CompScienceMeshes


Γ = CompScienceMeshes.meshrectangle(1.0, 1.0, 0.5)

X1 = BEAST.raviartthomas(Γ)

X = X1 × X1
ax = BEAST.NestedUnitRanges.nestedrange(X, 1, numfunctions)

coeffs = rand(numfunctions(X))
coeffs = BEAST.BlockArrays.BlockedVector(coeffs, (ax,))

u = BEAST.FEMFunction(coeffs, X)

@hilbertspace m[1:2]

u1 = u[m[1]]
u2 = u[m[2]]
u2 = -u1

@test u1 isa BEAST.FEMFunction
@test u2 isa BEAST.FEMFunction

@test u1.coeffs == -u2.coeffs