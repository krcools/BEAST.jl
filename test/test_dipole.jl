using Test

using CompScienceMeshes
using BEAST
using StaticArrays
using LinearAlgebra
using BlockArrays

for U in [Float32,Float64]

    c = U(3e8)
    μ0 = U(4*π*1e-7)
    μr = U(1.0)
    μ = μ0*μr
    ε0 = U(8.854187812e-12)
    εr = U(5)
    ε = ε0*εr
    c = U(1)/sqrt(ε*μ)
    local f = U(5e7)
    λ = c/f
    k = U(2*π/λ)
    local ω = k*c
    η = sqrt(μ/ε)

    a = U(1)
    Γ_orig = CompScienceMeshes.meshcuboid(a,a,a,U(0.1))
    local Γ = translate(Γ_orig,SVector(U(-a/2),U(-a/2),U(-a/2)))

    Φ, Θ = U.([0.0]), range(U(0),stop=U(π),length=100)
    pts = [point(U,cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]

    # This is an electric dipole
    # The pre-factor (1/ε) is used to resemble
    # (9.18) in Jackson's Classical Electrodynamics
    local E = U(1/ε) * dipolemw3d(location=SVector(U(0.4),U(0.2),U(0)),
                        orientation=U(1e-9).*SVector(U(0.5),U(0.5),U(0)),
                        wavenumber=k)

    local n = BEAST.NormalVector()

    𝒆 = (n × E) × n
    local H = (-1/(im*μ*ω))*curl(E)
    𝒉 = (n × H) × n

    𝓣 = Maxwell3D.singlelayer(wavenumber=k, alpha=-im*ω*μ, beta=-1/(im*ω*ε))
    𝓝 = BEAST.NCross()
    𝓚 = Maxwell3D.doublelayer(wavenumber=k)

    local X = raviartthomas(Γ)
    local Y = buffachristiansen(Γ)

    local T = Matrix(assemble(𝓣,X,X))
    e = Vector(assemble(𝒆,X))
    j_EFIE = T\e

    nf_E_EFIE = potential(MWSingleLayerField3D(𝓣), pts, j_EFIE, X)
    nf_H_EFIE = potential(BEAST.MWDoubleLayerField3D(𝓚), pts, j_EFIE, X)
    ff_E_EFIE = potential(MWFarField3D(𝓣), pts, j_EFIE, X)
    ff_H_EFIE = potential(BEAST.MWDoubleLayerFarField3D(𝓚), pts, j_EFIE, X)
    ff_H_EFIE_rotated = potential(n × BEAST.MWDoubleLayerFarField3D(𝓚), pts, -j_EFIE, n × X)
    ff_H_EFIE_doublerotated = potential(n × BEAST.MWDoubleLayerRotatedFarField3D(n × 𝓚), pts, -j_EFIE, X)


    @test norm(nf_E_EFIE - E.(pts))/norm(E.(pts)) ≈ 0 atol=0.01
    @test norm(nf_H_EFIE - H.(pts))/norm(H.(pts)) ≈ 0 atol=0.01
    @test norm(ff_E_EFIE - E.(pts, isfarfield=true))/norm(E.(pts, isfarfield=true)) ≈ 0 atol=0.001
    @test norm(ff_H_EFIE - H.(pts, isfarfield=true))/norm(H.(pts, isfarfield=true)) ≈ 0 atol=0.001
    @test norm(ff_H_EFIE_rotated - H.(pts, isfarfield=true))/norm(H.(pts, isfarfield=true)) ≈ 0 atol=0.001
    @test norm(ff_H_EFIE_doublerotated - H.(pts, isfarfield=true))/norm(H.(pts, isfarfield=true)) ≈ 0 atol=0.001
    @test ff_H_EFIE ≈ ff_H_EFIE_rotated rtol=1e-7
    @test ff_H_EFIE_rotated ≈ ff_H_EFIE_doublerotated rtol=1e-7


    K_bc = Matrix(assemble(𝓚,Y,X))
    G_nxbc_rt = Matrix(assemble(𝓝,Y,X))
    h_bc = Vector(assemble(𝒉,Y))
    M_bc = -U(0.5)*G_nxbc_rt + K_bc
    j_BCMFIE = M_bc\h_bc

    nf_E_BCMFIE = potential(MWSingleLayerField3D(𝓣), pts, j_BCMFIE, X)
    nf_H_BCMFIE = potential(BEAST.MWDoubleLayerField3D(𝓚), pts, j_BCMFIE, X)
    ff_E_BCMFIE = potential(MWFarField3D(𝓣), pts, j_BCMFIE, X)

    @test norm(nf_E_BCMFIE - E.(pts))/norm(E.(pts)) ≈ 0 atol=0.01
    @test norm(nf_H_BCMFIE - H.(pts))/norm(H.(pts)) ≈ 0 atol=0.01
    @test norm(ff_E_BCMFIE - E.(pts, isfarfield=true))/norm(E.(pts, isfarfield=true)) ≈ 0 atol=0.01

    H = dipolemw3d(location=SVector(U(0.0),U(0.0),U(0.3)),
                orientation=U(1e-9).*SVector(U(0.5),U(0.5),U(0)),
                wavenumber=k)

    # This time, we do not specify alpha and beta
    # We include η in the magnetic RHS
    𝓣 = Maxwell3D.singlelayer(wavenumber=k)

    𝒉 = (n × H) × n
    E = (1/(im*ε*ω))*curl(H)
    𝒆 = (n × E) × n

    X = raviartthomas(Γ)
    Y = buffachristiansen(Γ)

    T = Matrix(assemble(𝓣,X,X))
    e = Vector(assemble(𝒆,X))
    j_EFIE = T\e

    nf_E_EFIE = potential(MWSingleLayerField3D(wavenumber=k), pts, j_EFIE, X)
    nf_H_EFIE = potential(BEAST.MWDoubleLayerField3D(wavenumber=k), pts, j_EFIE, X) ./ η
    ff_E_EFIE = potential(MWFarField3D(wavenumber=k), pts, j_EFIE, X)

    @test norm(nf_E_EFIE - E.(pts))/norm(E.(pts)) ≈ 0 atol=0.01
    @test norm(nf_H_EFIE - H.(pts))/norm(H.(pts)) ≈ 0 atol=0.01
    @test norm(ff_E_EFIE - E.(pts, isfarfield=true))/norm(E.(pts, isfarfield=true)) ≈ 0 atol=0.01

    @hilbertspace s
    @hilbertspace t

    MFIE_discretebilform = @discretise(
        -U(0.5)*𝓝[t, s] + 𝓚[t, s] == η*𝒉[t],
        s∈X, t∈Y
    )

    j_BCMFIE = solve(MFIE_discretebilform)[Block(1)]

    nf_E_BCMFIE = potential(MWSingleLayerField3D(wavenumber=k), pts, j_BCMFIE, X)
    nf_H_BCMFIE = potential(BEAST.MWDoubleLayerField3D(wavenumber=k), pts, j_BCMFIE, X) ./ η
    ff_E_BCMFIE = potential(MWFarField3D(wavenumber=k), pts, j_BCMFIE, X)

    @test norm(j_BCMFIE - j_EFIE)/norm(j_EFIE) ≈ 0 atol=0.02
    @test norm(nf_E_BCMFIE - E.(pts))/norm(E.(pts)) ≈ 0 atol=0.01
    @test norm(nf_H_BCMFIE - H.(pts))/norm(H.(pts)) ≈ 0 atol=0.01
    @test norm(ff_E_BCMFIE - E.(pts, isfarfield=true))/norm(E.(pts, isfarfield=true)) ≈ 0 atol=0.01
end
