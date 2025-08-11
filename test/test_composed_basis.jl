######## test multiplied basis
using CompScienceMeshes
using LinearAlgebra
using BEAST
using Test
using StaticArrays

Γ = meshcuboid(1.0,1.0,1.0,1.0; generator=:compsciencemeshes)
X = raviartthomas(Γ)
func = BEAST.FunctionWrapper{Float64}(x->2.0)
Y = BEAST._BasisTimes(X,func)
K = Maxwell3D.doublelayer(wavenumber=1.0)
m1 = assemble(K,Y,Y)
m2 = assemble(K,X,X)
@test m1 ≈ 4*m2


L = lagrangecxd0(Γ)
Z = BEAST._BasisTimes(L,n)
m3 = assemble(K,Z,Z)
@test norm(m3) ≈ 0.10751130304639696

m1 = assemble(NCross(),X,X)
m2 = assemble(Identity(),X×n,X)

@test m1 ≈ m2
