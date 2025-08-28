using BEAST
using CompScienceMeshes
using StaticArrays
using Test

let
    ε0 = 8.854187821e-12
    μ0 = 4π*1e-7
    c0 = 1/sqrt(ε0*μ0)
    η0 = sqrt(μ0/ε0)

    f = 1e7 # 1 GHz
    ω = 2π*f
    λ = c0/f
    k = 2π/λ

    Γ = CompScienceMeshes.meshsegment(1.0,1.0)
    Γ = CompScienceMeshes.Mesh([SVector(0.0, 0.0), SVector(1.0, 1.0)], Γ.faces)

    X0 = lagrangecxd0(Γ)

    E0 = 1.0 # amplitude
    Einc = Helmholtz2D.planewave(; amplitude=E0, wavenumber=k, direction=SVector(1.0, 0.0))

    Hinc = curl(Einc)
    #hTM = assemble(TangentTrace(Hinc), X0)

    h1 = assemble(TangentTrace(Hinc), X0)
    h2 = assemble(BEAST.NormalDerivative(Einc), X0)

    @test norm(h1-h2) /norm(h2) ≈ 0  atol=1e-15

    ##
    Γ = CompScienceMeshes.meshsegment(1.0,0.5)

    X1 = lagrangec0d1(Γ)

    E0 = 1.0 # amplitude
    Einc = Helmholtz2D.planewave(; amplitude=E0, wavenumber=k, direction=SVector(0.0, 1.0))
    Hinc = curl(Einc)

    h1 = assemble(TangentTrace(Hinc), X1)
    h2 = assemble(BEAST.NormalDerivative(Einc), X1)

    @test norm(h1-h2) /norm(h2) ≈ 0
end