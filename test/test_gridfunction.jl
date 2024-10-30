@testitem "gridfunction" begin

# using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
# using Test
U = Float64

a = U(1)

# Cube between (3,3,3) and (4,4,4)
Γ = translate(CompScienceMeshes.meshcuboid(a,a,a,U(1.0)), SVector(3.0,3.0,3.0))

C0 = lagrangec0d1(Γ)

# We interpolate the function f(x, y, z) = z
coeffs = [v[3] for v in vertices(Γ)]

gfc0 = BEAST.FEMFunction(coeffs, C0)

# We integrate on a cube between the vertices (3,3,3) and (4,4,4)
# the  function f(x, y, z) = z
# L1 result is = 4 + 3 + 4 * (16/2 - 9/2) = 21
@test BEAST.Lp_integrate(gfc0; p=1) ≈ 21

# L2 result is sqrt(16 + 9 + 4 * (4^3/3 - 3^3/2))
@test BEAST.Lp_integrate(gfc0) ≈ sqrt(16 + 9 + 4 * (4^3/3 - 3^3/3))

#lincombgfs = BEAST.LinearCombinationOfAbstractMeshFunctions([1.0; -0.5], [gfc0, gfc0])
lincombgfs = -0.5*gfc0 + gfc0

# We integrate on a cube between the vertices (3,3,3) and (4,4,4)
# the  function f(x, y, z) = z
# L1 result is = 4 + 3 + 4 * (16/2 - 9/2) = 21
@test BEAST.Lp_integrate(lincombgfs; p=1) ≈ 21/2

# L2 result is sqrt(16 + 9 + 4 * (4^3/3 - 3^3/2))
@test BEAST.Lp_integrate(lincombgfs) ≈ sqrt(16 + 9 + 4 * (4^3/3 - 3^3/3))/2

CX = lagrangecxd0(Γ)

# We interpolate the function f(x, y, z) = 2
# with piecewise constant lagrange elements
constgfcx = BEAST.FEMFunction([2.0 for c in cells(Γ)], CX)
@test BEAST.Lp_integrate(constgfcx; p=1) ≈ 12
@test BEAST.Lp_integrate(constgfcx) ≈ sqrt(24)

# We interpolate the function f(x, y, z) = z
# with piecewise constant lagrange elements
coeffsCx = [((Γ.vertices[c[1]] + Γ.vertices[c[2]] + Γ.vertices[c[3]])/3)[3] for c in cells(Γ)]

gfcx = BEAST.FEMFunction(coeffsCx, CX)

# Even though the function, should be able to interpolate it
@test BEAST.Lp_integrate(gfcx; p=1) ≈ 21

# Compute the difference in approximation between constant and linear elements
#difflincomb1 = BEAST.LinearCombinationOfAbstractMeshFunctions([1.0; -1.0], [gfcx, gfc0])
difflincomb1 = gfcx - gfc0

@test BEAST.Lp_integrate(difflincomb1)/BEAST.Lp_integrate(gfc0) ≈ 0.03866223365135301

## Test GlobalFunction
f(x) = x[3]
glbf = BEAST.GlobalFunction(f, Γ, Vector(1:numcells(Γ)))

@test BEAST.Lp_integrate(glbf; p=1) ≈ 21

g(x) = x[3]^2
glbf2 = BEAST.GlobalFunction(g, Γ, Vector(1:numcells(Γ)))

@test BEAST.Lp_integrate(glbf2; p=1) ≈ 74.33333333333334

## Test linear combination of global and fem functions
idxfem, idxglobal = BEAST.indices_splitfemglobal(gfc0 - glbf + 0.0gfcx)

@test length(idxfem) == 2
@test length(idxglobal) == 1

##
@test BEAST.Lp_integrate(gfc0 - glbf) / BEAST.Lp_integrate(glbf) ≈ 0 atol=1e-16

@test BEAST.Lp_integrate(gfcx - glbf) / BEAST.Lp_integrate(glbf) ≈ 0.038662233651353003

## Global function on limited supported
f(x) = x[3]

idx_topplate = Int64[]

for (i, face) in enumerate(cells(Γ))
    center = (Γ.vertices[face[1]] + Γ.vertices[face[2]] + Γ.vertices[face[3]]) / 3.0
    if center[3] > 3.99
        push!(idx_topplate, i)
    end
end

glbf_toppplate = BEAST.GlobalFunction(f, Γ, idx_topplate)
@test BEAST.Lp_integrate(glbf_toppplate, p=1) ≈ 4.000000000000003

@test BEAST.Lp_integrate(gfc0 - glbf_toppplate, p=1) ≈ 21 - 4 # Topplate accounts for 4

##

Γ2 = CompScienceMeshes.meshcuboid(a,a,a,U(1.0))

coeffsΓ2 = [v[3] for v in vertices(Γ2)]

C0Γ2 = lagrangec0d1(Γ2)

gfc0Γ2 = BEAST.FEMFunction(coeffsΓ2, C0Γ2)

@test_throws ErrorException("Functions must be defined over the same geometry.") gfc0 + gfc0Γ2

@test_throws ErrorException("Functions must be defined over the same geometry.") lincombgfs + 1.0gfc0Γ2

## Test mixed combinations
@test BEAST.Lp_integrate(gfc0 + 1im * gfc0; p=1) ≈ sqrt(2) * 21

end