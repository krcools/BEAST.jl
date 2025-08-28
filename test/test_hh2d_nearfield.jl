using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test

@testset "2D-Helmholtz potential operators" begin
    # r = 10.0
    # λ = 20 * r
    # k = 2 * π / λ
    # sphere = meshsphere(r, 0.2 * r)
    k = 0.031415926535897934
    r = 10.0
    circle = CompScienceMeshes.meshcircle(r, 2.0)

    X0 = lagrangecxd0(circle)
    X1 = lagrangec0d1(circle)

    S = Helmholtz2D.singlelayer(; gamma=im * k)
    D = Helmholtz2D.doublelayer(; gamma=im * k)
    Dt = Helmholtz2D.doublelayer_transposed(; gamma=im * k)
    N = Helmholtz2D.hypersingular(; gamma=im * k)

    q = 100.0
    ϵ = 1.0

    # Interior problem
    # Formulations from Sauter and Schwab, Boundary Element Methods(2011), Chapter 3.4.1.1
  
    pos1 = SVector(r * 1.5, 0.0)  # positioning of point charges
    pos2 = SVector(-r * 1.5, 0.0)

    charge1 = Helmholtz2D.monopole(position=pos1, amplitude=q/(4*π*ϵ), wavenumber=k)
    charge2 = Helmholtz2D.monopole(position=pos2, amplitude=-q/(4*π*ϵ), wavenumber=k)

    # Potential of point charges

    Φ_inc(x) = charge1(x) + charge2(x)

    gD0 = assemble(DirichletTrace(charge1), X0) + assemble(DirichletTrace(charge2), X0)
    gD1 = assemble(DirichletTrace(charge1), X1) + assemble(DirichletTrace(charge2), X1)
    gN = assemble(∂n(charge1), X1) + assemble(BEAST.n ⋅ grad(charge2), X1)

    G = assemble(Identity(), X1, X1)
    o = ones(numfunctions(X1))

    # Interior Dirichlet problem - compare Sauter & Schwab eqs. 3.81
    M_IDPSL = assemble(S, X0, X0) # Single layer (SL)
    M_IDPDL = (-1 / 2 * assemble(Identity(), X1, X1) + assemble(D, X1, X1)) # Double layer (DL)

    # Interior Neumann problem
    # Neumann derivative from DL potential with deflected nullspace
    M_INPDL = assemble(N, X1, X1)
    # Neumann derivative from SL potential with deflected nullspace
    M_INPSL = (1 / 2 * assemble(Identity(), X1, X1) + assemble(Dt, X1, X1))

    ρ_IDPSL = M_IDPSL \ (-gD0)
    ρ_IDPDL = M_IDPDL \ (-gD1)

    ρ_INPDL = M_INPDL \ (gN)
    ρ_INPSL = M_INPSL \ (-gN)

    pts = meshcircle(0.8 * r, 0.8 * 0.6 * r).vertices # sphere inside on which the potential and field are evaluated

    pot_IDPSL = potential(HH2DSingleLayerNear(S), pts, ρ_IDPSL, X0; type=ComplexF64)
    pot_IDPDL = potential(HH2DDoubleLayerNear(D), pts, ρ_IDPDL, X1; type=ComplexF64)

    pot_INPSL = potential(HH2DSingleLayerNear(S), pts, ρ_INPSL, X1; type=ComplexF64)
    pot_INPDL = potential(HH2DDoubleLayerNear(D), pts, ρ_INPDL, X1; type=ComplexF64)

    # Total field inside should be zero
    err_IDPSL_pot = norm(pot_IDPSL + Φ_inc.(pts)) / norm(Φ_inc.(pts))
    err_IDPDL_pot = norm(pot_IDPDL + Φ_inc.(pts)) / norm(Φ_inc.(pts))
    err_INPSL_pot = norm(pot_INPSL + Φ_inc.(pts)) / norm(Φ_inc.(pts))
    err_INPDL_pot = norm(pot_INPDL + Φ_inc.(pts)) / norm(Φ_inc.(pts))

    # Efield(x) = - grad Φ_inc(x)
    Efield(x) =  -grad(charge1)(x) + -grad(charge2)(x)

    field_IDPSL = -potential(HH2DDoubleLayerTransposedNear(Dt), pts, ρ_IDPSL, X0; type=SVector{2,ComplexF64})
    field_IDPDL = -potential(HH2DHyperSingularNear(N), pts, ρ_IDPDL, X1; type=SVector{2,ComplexF64})
    field_INPSL = -potential(HH2DDoubleLayerTransposedNear(Dt), pts, ρ_INPSL, X1; type=SVector{2,ComplexF64})
    field_INPDL = -potential(HH2DHyperSingularNear(N), pts, ρ_INPDL, X1; type=SVector{2,ComplexF64})

    err_IDPSL_field = norm(field_IDPSL + Efield.(pts)) / norm(Efield.(pts))
    err_IDPDL_field = norm(field_IDPDL + Efield.(pts)) / norm(Efield.(pts))
    err_INPSL_field = norm(field_INPSL + Efield.(pts)) / norm(Efield.(pts))
    err_INPDL_field = norm(field_INPDL + Efield.(pts)) / norm(Efield.(pts))

    # Exterior problem 
    # formulations from Sauter and Schwab, Boundary Element Methods(2011), Chapter 3.4.1.2

    pos1 = SVector(r * 0.5, 0.0)
    pos2 = SVector(-r * 0.5, 0.0)

    charge1 = Helmholtz2D.monopole(position=pos1, amplitude=q/(4*π*ϵ), wavenumber=k)
    charge2 = Helmholtz2D.monopole(position=pos2, amplitude=-q/(4*π*ϵ), wavenumber=k)

    gD0 = assemble(DirichletTrace(charge1), X0) + assemble(DirichletTrace(charge2), X0)
    gD1 = assemble(DirichletTrace(charge1), X1) + assemble(DirichletTrace(charge2), X1)
    gN = assemble(∂n(charge1), X1) + assemble(∂n(charge2), X1)

    G = assemble(Identity(), X1, X1)
    o = ones(numfunctions(X1))

    M_EDPSL = assemble(S, X0, X0)
    M_EDPDL = (1 / 2 * assemble(Identity(), X1, X1) + assemble(D, X1, X1))

    M_ENPDL = assemble(N, X1, X1)
    M_ENPSL = -1 / 2 * assemble(Identity(), X1, X1) + assemble(Dt, X1, X1)

    ρ_EDPSL = M_EDPSL \ (-gD0)
    ρ_EDPDL = M_EDPDL \ (-gD1)

    ρ_ENPDL = M_ENPDL \ gN
    ρ_ENPSL = M_ENPSL \ (-gN)

    testcircle = meshcircle(1.2 * r, 1.2 * 0.6 * r)
    pts = testcircle.vertices[norm.(testcircle.vertices) .> r]

    pot_EDPSL = potential(HH2DSingleLayerNear(S), pts, ρ_EDPSL, X0; type=ComplexF64)
    pot_EDPDL = potential(HH2DDoubleLayerNear(D), pts, ρ_EDPDL, X1; type=ComplexF64)
    pot_ENPDL = potential(HH2DDoubleLayerNear(D), pts, ρ_ENPDL, X1; type=ComplexF64)
    pot_ENPSL = potential(HH2DSingleLayerNear(S), pts, ρ_ENPSL, X1; type=ComplexF64)

    err_EDPSL_pot = norm(pot_EDPSL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))
    err_EDPDL_pot = norm(pot_EDPDL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))
    err_ENPSL_pot = norm(pot_ENPSL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))
    err_ENPDL_pot = norm(pot_ENPDL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))

    field_EDPSL = -potential(HH2DDoubleLayerTransposedNear(Dt), pts, ρ_EDPSL, X0; type=SVector{2,ComplexF64})
    field_EDPDL = -potential(HH2DHyperSingularNear(N), pts, ρ_EDPDL, X1; type=SVector{2,ComplexF64})
    field_ENPSL = -potential(HH2DDoubleLayerTransposedNear(Dt), pts, ρ_ENPSL, X1; type=SVector{2,ComplexF64})
    field_ENPDL = -potential(HH2DHyperSingularNear(N), pts, ρ_ENPDL, X1; type=SVector{2,ComplexF64})

    err_EDPSL_field = norm(field_EDPSL + Efield.(pts)) / norm(Efield.(pts))
    err_EDPDL_field = norm(field_EDPDL + Efield.(pts)) / norm(Efield.(pts))
    err_ENPSL_field = norm(field_ENPSL + Efield.(pts)) / norm(Efield.(pts))
    err_ENPDL_field = norm(field_ENPDL + Efield.(pts)) / norm(Efield.(pts))

    # errors of interior problems
    @test err_IDPSL_pot < 1e-3
    @test err_IDPDL_pot < 1e-3
    @test err_INPSL_pot < 1e-3
    @test err_INPDL_pot < 1e-3

    @test err_IDPSL_field < 1e-2
    @test err_IDPDL_field < 1e-2
    @test err_INPSL_field < 1e-2
    @test err_INPDL_field < 1e-2

    # errors of exterior problems
    @test err_EDPSL_pot < 1e-3
    @test err_EDPDL_pot < 1e-3
    @test err_ENPSL_pot < 1e-3
    @test err_ENPDL_pot < 1e-3

    @test err_EDPSL_field < 1e-2
    @test err_EDPDL_field < 1e-2
    @test err_ENPSL_field < 1e-2
    @test err_ENPDL_field < 1e-2
end