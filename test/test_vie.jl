using BEAST
using CompScienceMeshes
using LinearAlgebra

using Test


"""
Tests the structure of the VIE implementation (function spaces, numerical integration, traces, ...) 
by using the divergence theorem: ∫∫∫_Ω ∫∫∫_Ω ... = ∫∫∫_Ω ∫∫∫_Ω ... + ∫∫_Γ ∫∫∫_Ω ...
This test does not include the Operators of the EVIE, DVIE and one operator of the LSVIE.
"""


@testset "Volume Integral Equations" begin 

    qs = BEAST.SauterSchwab3DQStrat(3,3,3,3,3,3)

    r = 1.0
    h = 0.4
    mesh = CompScienceMeshes.tetmeshsphere(r,h)
    bnd = boundary(mesh)

    function tau(x)
        return 1.0 + 2.0*x[1]^2*(x[2] - 1.0) + (x[3] - 0.5)^2 - 7.0*x[1]*x[2]*x[3]
    end

    #using PlotlyJS
    #PlotlyJS.plot(CompScienceMeshes.wireframe(mesh))

    ## LSVIE
    X = lagrangec0d1(mesh; dirichlet = false)
    strc = X -> strace(X,bnd)

    V, B = VIE.hhvolume(tau = tau, wavenumber = 1.0), VIE.hhboundary(tau = tau, wavenumber = 1.0)
    Y = VIE.hhvolumegradG(tau = tau, wavenumber = 1.0)

    M_a = assemble(Y, X, X, quadstrat = qs)
    M_b = assemble(B, strc(X), X, quadstrat = qs) + assemble(V, X, X, quadstrat = qs)
    err = norm(M_a - M_b)/norm(M_b)
    @show err
    @test err < 0.005


    ## DVIE
    X = nedelecd3d(mesh)
    ntrc = X -> ntrace(X,bnd)
    
    TestOp1 = BEAST.DVIE_TestOp1(1.0, 1.0, tau)
    TestOp2 = BEAST.DVIE_TestOp2(1.0, 1.0, tau)
    TestOp3 = BEAST.DVIE_TestOp3(1.0, -1.0, tau)

    M_a = assemble(TestOp1, X, X, quadstrat = qs)
    M_b = assemble(TestOp2, ntrc(X), X, quadstrat = qs) + assemble(TestOp3, X, X, quadstrat = qs)
    
    err = norm(M_a - M_b)/norm(M_b)
    @show err
    @test err < 0.02


    ## EVIE
    X = nedelecc3d(mesh)
    ttrc = X -> ttrace(X,bnd)

    TestOp1 = BEAST.EVIE_TestOp1(1.0, 1.0, tau)
    TestOp2 = BEAST.EVIE_TestOp2(1.0, 1.0, tau)
    TestOp3 = BEAST.EVIE_TestOp3(1.0, 1.0, tau)

    M_a = assemble(TestOp1, X, X, quadstrat = qs)
    M_b = assemble(TestOp2, ttrc(X), X, quadstrat = qs) + assemble(TestOp3, X, X, quadstrat = qs)
    
    err = norm(M_a - M_b)/norm(M_b)
    @show err
    @test err < 0.01

end