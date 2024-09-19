######## test multiplied basis
using CompScienceMeshes
using LinearAlgebra
using BEAST
using Test

Γ = meshcuboid(1.0,1.0,1.0,1.0)
X = raviartthomas(Γ)
f(x) = 2.0
Y = BEAST._BasisTimes(X,f)
K = Maxwell3D.doublelayer(wavenumber=1.0)
m1 = assemble(K,Y,Y)
m2 = assemble(K,X,X)
@test m1 ≈ 4*m2

U = BEAST._BasisTimes(BEAST._BasisDot(x->(@SVector [1.0,1.0,1.0]),n),X)
m = assemble(K,U,U)

L = lagrangecxd0(Γ)
Z = BEAST._BasisTimes(L,n)
m3 = assemble(K,Z,Z)
@test norm(m3) ≈ 0.08420116178577139



##### go in code 
# using StaticArrays
# s = simplex((@SVector [1.0,0.0,0.0]),(@SVector [0.0,1.0,0.0]),(@SVector [0.0,0.0,0.0]))
# p = neighborhood(s,[0.5,0.2])

# rtref = BEAST.RTRefSpace{Float64}()
# f(x)= 2.0

# mref = BEAST._LocalBasisTimes(Float64,BEAST.FunctionWrapper(f),rtref)
