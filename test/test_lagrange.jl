using Test
using BEAST
using CompScienceMeshes

@testset "Higher order 1D Lagrange Elements" begin
    line = meshsegment(2.0, 0.5)
    numsegments = CompScienceMeshes.numcells(line)
    g = BEAST.ScalarTrace{Float64}(x -> x[1]^5)
    glbf = BEAST.GlobalFunction(g, line, Vector(1:numcells(line)))
    @test BEAST.Lp_integrate(glbf; p=1) ≈ 10.66666666666


    @testset "C0" begin

        for order in 1:10
            C0 = lagrangec0(line, order=order, dirichlet=false)
            @test numfunctions(C0) == (order-1)*numsegments + numsegments + 1
        end

        # We integrate x^5 from 0 to 2
        # Analytic outcome: 32/3

        # Order 1
        C0 = lagrangec0(line, order=1, dirichlet=false)
        coeffs = DofInterpolate(C0, g) 
        gfc0 = BEAST.FEMFunction(coeffs, C0)
        @test BEAST.Lp_integrate(gfc0; p=1) ≈ 12.3125

        # Order 2
        C0 = lagrangec0(line, order=2, dirichlet=false)
        coeffs = DofInterpolate(C0, g) 
        gfc0 = BEAST.FEMFunction(coeffs, C0)
        @test BEAST.Lp_integrate(gfc0; p=1) ≈ 10.672838520122093

        # Order 3
        C0 = lagrangec0(line, order=3, dirichlet=false)
        coeffs = DofInterpolate(C0, g) 
        gfc0 = BEAST.FEMFunction(coeffs, C0)
        @test BEAST.Lp_integrate(gfc0; p=1) ≈ 10.6689814814888

        @test C0.pos[6][1] ≈ 0.166666666666666

        # Order 4
        C0 = lagrangec0(line, order=4, dirichlet=false)
        coeffs = DofInterpolate(C0, g) 
        gfc0 = BEAST.FEMFunction(coeffs, C0)
        @test BEAST.Lp_integrate(gfc0; p=1) ≈ 10.66668405

        # Order 5
        C0 = lagrangec0(line, order=5, dirichlet=false)
        coeffs = DofInterpolate(C0, g) 
        gfc0 = BEAST.FEMFunction(coeffs, C0)
        @test BEAST.Lp_integrate(gfc0; p=1) ≈ 10.6666666666

        @test C0.pos[6][1] ≈ 0.1
    end
    ##

    @testset "CX" begin

        for order in 1:10
            CX = lagrangecx(line, order=order)
            @test numfunctions(CX) == (order+1)*numsegments 
        end

        # Order 1
        CX = lagrangecx(line, order=1)
        coeffs = DofInterpolate(CX, g) 
        gfc0 = BEAST.FEMFunction(coeffs, CX)
        @test BEAST.Lp_integrate(gfc0; p=1) ≈ 12.3125

        CX = lagrangecx(line, order=2)
        coeffs = DofInterpolate(CX, g) 
        gfc0 = BEAST.FEMFunction(coeffs, CX)
        @test BEAST.Lp_integrate(gfc0; p=1) ≈ 10.672838520122093

        CX = lagrangecx(line, order=5)
        coeffs = DofInterpolate(CX, g) 
        gfc0 = BEAST.FEMFunction(coeffs, CX)
        @test BEAST.Lp_integrate(gfc0; p=1) ≈ 10.6666666666

        @test CX.pos[9] ≈ point(0.1, 0.0)
    end
end