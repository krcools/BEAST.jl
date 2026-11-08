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

            # Which kernel values are evaluated: the boundary integral operators evaluate
            # what their integrand reads, the near-field operators both.
            for (op, evaluatesgreen, evaluatesgradgreen) in (
                (BEAST.HH2DSingleLayerFDBIO(1.0, γ), true, false),
                (BEAST.HH2DDoubleLayerFDBIO(1.0, γ), false, true),
                (BEAST.HH2DDoubleLayerTransposedFDBIO(1.0, γ), false, true),
                (BEAST.HH2DHyperSingularFDBIO(1.0, 1.0, γ), true, false),
                (BEAST.HH2DSingleLayerNear(1.0, γ), true, true),
                (BEAST.HH2DDoubleLayerNear(1.0, γ), true, true),
                (BEAST.HH2DDoubleLayerTransposedNear(1.0, γ), true, true),
                (BEAST.HH2DHyperSingularNear(1.0, 1.0, γ), true, true),
            )
                # Type-stable: the traits that select the kernel values are constant
                kv = @inferred kernelvals(op, mp1, mp2)
                @test kv.vect ≈ r
                @test kv.dist ≈ R

                # A kernel value that is not evaluated is `nothing`: the Hankel or Bessel
                # function is computed only when it is read.
                if evaluatesgreen
                    @test kv.green ≈ green
                else
                    @test isnothing(kv.green)
                end
                if evaluatesgradgreen
                    @test kv.gradgreen ≈ gradgreen
                else
                    @test isnothing(kv.gradgreen)
                end

                # The kernel values are evaluated at every quadrature point of an
                # assembly and must not allocate.
                @test isbits(kv)
                allocations(op, mp1, mp2)
                @test allocations(op, mp1, mp2) <= specialallocations
            end
        end
    end

    @testset "assembly reads only the evaluated kernel values" begin
        # Coincident, adjacent and distant segments of a closed curve, including the
        # singular quadrature rules: an integrand that read a kernel value that is not
        # evaluated (`nothing`) would throw. In single precision, the Hankel and Bessel
        # functions return double precision.
        for F in (Float64, Float32)
            X = lagrangec0d1(meshcircle(F(1.0), F(0.5)))
            for γ in (Val(0), F(1.5), F(1.5) * im, complex(F(1.0), F(0.5)))
                for op in (
                    BEAST.HH2DSingleLayerFDBIO(F(1.0), γ),
                    BEAST.HH2DDoubleLayerFDBIO(F(1.0), γ),
                    BEAST.HH2DDoubleLayerTransposedFDBIO(F(1.0), γ),
                    BEAST.HH2DHyperSingularFDBIO(F(1.0), F(1.0), γ),
                )
                    @test all(isfinite, assemble(op, X, X))
                end
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