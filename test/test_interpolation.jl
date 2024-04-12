using BEAST
using CompScienceMeshes
using Test

Γ =  meshrectangle(1.0, 1.0, 0.4)

XN = lagrangecxd0(Γ) # lagrangec0d1(Γₑ₁; dirichlet=true) # Dirichlet=true -> no boundary function
XD = lagrangec0d1(Γ)

mp = Helmholtz3D.monopole(position=SVector(0.0, 0.0, 3.0))
gradmp = grad(mp)

φ = DofInterpolate(XD, mp)

@test eltype(φ) == Float64

for i in eachindex(φ)
    @test mp(XD.pos[i]) == φ[i] 
end

σ = DofInterpolate(XN, gradmp)

@test eltype(σ) == Float64

for i in eachindex(σ)
    # Orientation of meshrectangle is in -z direction
    @test dot(SVector(0.0, 0.0, -1.0), gradmp(XN.pos[i])) == σ[i]
end