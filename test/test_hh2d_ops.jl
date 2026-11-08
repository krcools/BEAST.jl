using Test
using CompScienceMeshes
using BEAST
using LinearAlgebra
using SpecialFunctions
using StaticArrays

@testset "HelmholtzOperator2D" begin

    v1 = SVector(1.0, 0.0)
    v2 = SVector(0.0, 0.0)

    verts = [v1, v2]
    faces = [CompScienceMeshes.SimplexGraph(SVector(1,2))]

    m = CompScienceMeshes.Mesh(verts, faces)

    ch = chart(m, 1)

    mp1 = neighborhood(ch, SVector(0.15))
    mp2 = neighborhood(ch, SVector(0.85))

    @testset "hh2d_makegammacomplexifneeded" begin
        # Test positive real gamma
        γ_pos = 1.0
        result = BEAST.hh2d_makegammacomplexifneeded(γ_pos)
        @test result == 1.0
        
        # Test negative real gamma (should be complexified)
        γ_neg = -1.0
        result = BEAST.hh2d_makegammacomplexifneeded(γ_neg)
        @test isa(result, Complex)
        @test real(result) == -1.0
        
        # Test complex gamma (should pass through)
        γ_complex = 1.0 + 2.0im
        result = BEAST.hh2d_makegammacomplexifneeded(γ_complex)
        @test result == γ_complex
    end
    
    @testset "kernelvals function test" begin
        α = 1.0
        γ = 1.0 + 1.0im
        op = BEAST.HH2DSingleLayerFDBIO(α, γ)
        @test op.alpha == α
        @test op.gamma == γ

        kv = kernelvals(op, mp1, mp2)

        @test kv.gamma == γ
    end

    @testset "kernelvals: values, no allocations" begin
        r = cartesian(mp1) - cartesian(mp2)
        R = norm(r)
        allocations(op, x, y) = @allocated kernelvals(op, x, y)
        # On Julia 1.10, besselk and hankelh2 allocate themselves (16 B per call for most
        # arguments); the kernel values must not allocate anything beyond that.
        besselkallocations(x) = @allocated (besselk(0, x), besselk(1, x))
        hankelallocations(x) = @allocated (hankelh2(0, x), hankelh2(1, x))

        # Laplace, modified Helmholtz, and Helmholtz with real and complex wavenumber
        for γ in (Val(0), 1.5, 1.5im, 1.0 + 0.5im)
            if γ isa Val
                green = -log(R) / 2π
                gradgreen = -r / (2π * R^2)
                specialallocations = 0
            elseif γ isa Real
                green = besselk(0, γ * R) / 2π
                gradgreen = -γ / 2π * besselk(1, γ * R) * r / R
                besselkallocations(γ * R)
                specialallocations = besselkallocations(γ * R)
            else
                k = iszero(real(γ)) ? imag(γ) : -im * γ
                green = -im / 4 * hankelh2(0, k * R)
                gradgreen = k * im / 4 * hankelh2(1, k * R) * r / R
                hankelallocations(k * R)
                specialallocations = hankelallocations(k * R)
            end

            for op in (
                BEAST.HH2DSingleLayerFDBIO(1.0, γ),
                BEAST.HH2DDoubleLayerFDBIO(1.0, γ),
                BEAST.HH2DDoubleLayerTransposedFDBIO(1.0, γ),
                BEAST.HH2DHyperSingularFDBIO(1.0, 1.0, γ),
                BEAST.HH2DSingleLayerNear(1.0, γ),
                BEAST.HH2DDoubleLayerNear(1.0, γ),
                BEAST.HH2DDoubleLayerTransposedNear(1.0, γ),
                BEAST.HH2DHyperSingularNear(1.0, 1.0, γ),
            )
                kv = kernelvals(op, mp1, mp2)
                @test kv.vect ≈ r
                @test kv.dist ≈ R
                @test kv.green ≈ green
                @test kv.gradgreen ≈ gradgreen

                # The kernel values are evaluated at every quadrature point of an
                # assembly and must not allocate.
                @test isbits(kv)
                allocations(op, mp1, mp2)
                @test allocations(op, mp1, mp2) <= specialallocations
            end
        end
    end
    
    @testset "scalartype" begin
        α = 1.0
        γ = 1.0 + 1.0im
        op = BEAST.HH2DSingleLayerFDBIO(α, γ)
        st = scalartype(op)
        @test st == promote_type(typeof(α), typeof(γ))
    end
    
    @testset "isstatic" begin
        α = 1.0
        γ = 1.0 + 1.0im
        op = BEAST.HH2DSingleLayerFDBIO(α, γ)
        @test !BEAST.isstatic(op)

        α = 1.0
        γ = Val(0)
        op = BEAST.HH2DSingleLayerFDBIO(α, γ)
        @test BEAST.isstatic(op)
    end
    
end