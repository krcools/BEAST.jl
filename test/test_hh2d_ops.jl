using Test
using CompScienceMeshes
using BEAST
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