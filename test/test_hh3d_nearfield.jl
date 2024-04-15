using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test

@testset "Helmholtz potential operators" begin
    # r = 10.0
    # λ = 20 * r
    # k = 2 * π / λ
    # sphere = meshsphere(r, 0.2 * r)
    k = 0.031415926535897934
    sphere = readmesh(joinpath(dirname(pathof(BEAST)),"../test/assets/sphere_rad=10_h=2.in"))

    X0 = lagrangecxd0(sphere)
    X1 = lagrangec0d1(sphere)

    S = Helmholtz3D.singlelayer(; gamma=im * k)
    D = Helmholtz3D.doublelayer(; gamma=im * k)
    Dt = Helmholtz3D.doublelayer_transposed(; gamma=im * k)
    N = Helmholtz3D.hypersingular(; gamma=im * k)

    q = 100.0
    ϵ = 1.0

    # Interior problem
    # Formulations from Sauter and Schwab, Boundary Element Methods(2011), Chapter 3.4.1.1
    r = 10.0
    pos1 = SVector(r * 1.5, 0.0, 0.0)  # positioning of point charges
    pos2 = SVector(-r * 1.5, 0.0, 0.0)

    charge1 = Helmholtz3D.monopole(position=pos1, amplitude=q/(4*π*ϵ), wavenumber=k)
    charge2 = Helmholtz3D.monopole(position=pos2, amplitude=-q/(4*π*ϵ), wavenumber=k)

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
    M_INPDL = assemble(N, X1, X1) + G * o * o' * G
    # Neumann derivative from SL potential with deflected nullspace
    M_INPSL = (1 / 2 * assemble(Identity(), X1, X1) + assemble(Dt, X1, X1)) + G * o * o' * G 

    ρ_IDPSL = M_IDPSL \ (-gD0)
    ρ_IDPDL = M_IDPDL \ (-gD1)

    ρ_INPDL = M_INPDL \ (gN)
    ρ_INPSL = M_INPSL \ (-gN)

    pts = meshsphere(0.8 * r, 0.8 * 0.6 * r).vertices # sphere inside on which the potential and field are evaluated

    pot_IDPSL = potential(HH3DSingleLayerNear(im * k), pts, ρ_IDPSL, X0; type=ComplexF64)
    pot_IDPDL = potential(HH3DDoubleLayerNear(im * k), pts, ρ_IDPDL, X1; type=ComplexF64)

    pot_INPSL = potential(HH3DSingleLayerNear(im * k), pts, ρ_INPSL, X1; type=ComplexF64)
    pot_INPDL = potential(HH3DDoubleLayerNear(im * k), pts, ρ_INPDL, X1; type=ComplexF64)

    # Total field inside should be zero
    err_IDPSL_pot = norm(pot_IDPSL + Φ_inc.(pts)) / norm(Φ_inc.(pts))
    err_IDPDL_pot = norm(pot_IDPDL + Φ_inc.(pts)) / norm(Φ_inc.(pts))
    err_INPSL_pot = norm(pot_INPSL + Φ_inc.(pts)) / norm(Φ_inc.(pts))
    err_INPDL_pot = norm(pot_INPDL + Φ_inc.(pts)) / norm(Φ_inc.(pts))

    # Efield(x) = - grad Φ_inc(x)
    Efield(x) =  -grad(charge1)(x) + -grad(charge2)(x)

    field_IDPSL = -potential(HH3DDoubleLayerTransposedNear(im * k), pts, ρ_IDPSL, X0)
    field_IDPDL = -potential(HH3DHyperSingularNear(im * k), pts, ρ_IDPDL, X1)
    field_INPSL = -potential(HH3DDoubleLayerTransposedNear(im * k), pts, ρ_INPSL, X1)
    field_INPDL = -potential(HH3DHyperSingularNear(im * k), pts, ρ_INPDL, X1)

    err_IDPSL_field = norm(field_IDPSL + Efield.(pts)) / norm(Efield.(pts))
    err_IDPDL_field = norm(field_IDPDL + Efield.(pts)) / norm(Efield.(pts))
    err_INPSL_field = norm(field_INPSL + Efield.(pts)) / norm(Efield.(pts))
    err_INPDL_field = norm(field_INPDL + Efield.(pts)) / norm(Efield.(pts))

    # Exterior problem 
    # formulations from Sauter and Schwab, Boundary Element Methods(2011), Chapter 3.4.1.2

    pos1 = SVector(r * 0.5, 0.0, 0.0)
    pos2 = SVector(-r * 0.5, 0.0, 0.0)

    charge1 = Helmholtz3D.monopole(position=pos1, amplitude=q/(4*π*ϵ), wavenumber=k)
    charge2 = Helmholtz3D.monopole(position=pos2, amplitude=-q/(4*π*ϵ), wavenumber=k)

    gD0 = assemble(DirichletTrace(charge1), X0) + assemble(DirichletTrace(charge2), X0)
    gD1 = assemble(DirichletTrace(charge1), X1) + assemble(DirichletTrace(charge2), X1)
    gN = assemble(∂n(charge1), X1) + assemble(∂n(charge2), X1)

    G = assemble(Identity(), X1, X1)
    o = ones(numfunctions(X1))

    M_EDPSL = assemble(S, X0, X0)
    M_EDPDL = (1 / 2 * assemble(Identity(), X1, X1) + assemble(D, X1, X1))

    M_ENPDL = assemble(N, X1, X1) + G * o * o' * G
    M_ENPSL = -1 / 2 * assemble(Identity(), X1, X1) + assemble(Dt, X1, X1) + G * o * o' * G

    ρ_EDPSL = M_EDPSL \ (-gD0)
    ρ_EDPDL = M_EDPDL \ (-gD1)

    ρ_ENPDL = M_ENPDL \ gN
    ρ_ENPSL = M_ENPSL \ (-gN)

    testsphere = meshsphere(1.2 * r, 1.2 * 0.6 * r)
    pts = testsphere.vertices[norm.(testsphere.vertices) .> r]

    pot_EDPSL = potential(HH3DSingleLayerNear(im * k), pts, ρ_EDPSL, X0; type=ComplexF64)
    pot_EDPDL = potential(HH3DDoubleLayerNear(im * k), pts, ρ_EDPDL, X1; type=ComplexF64)
    pot_ENPDL = potential(HH3DDoubleLayerNear(im * k), pts, ρ_ENPDL, X1; type=ComplexF64)
    pot_ENPSL = potential(HH3DSingleLayerNear(im * k), pts, ρ_ENPSL, X1; type=ComplexF64)

    err_EDPSL_pot = norm(pot_EDPSL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))
    err_EDPDL_pot = norm(pot_EDPDL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))
    err_ENPSL_pot = norm(pot_ENPSL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))
    err_ENPDL_pot = norm(pot_ENPDL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))

    field_EDPSL = -potential(HH3DDoubleLayerTransposedNear(im * k), pts, ρ_EDPSL, X0)
    field_EDPDL = -potential(HH3DHyperSingularNear(im * k), pts, ρ_EDPDL, X1)
    field_ENPSL = -potential(HH3DDoubleLayerTransposedNear(im * k), pts, ρ_ENPSL, X1)
    field_ENPDL = -potential(HH3DHyperSingularNear(im * k), pts, ρ_ENPDL, X1)

    err_EDPSL_field = norm(field_EDPSL + Efield.(pts)) / norm(Efield.(pts))
    err_EDPDL_field = norm(field_EDPDL + Efield.(pts)) / norm(Efield.(pts))
    err_ENPSL_field = norm(field_ENPSL + Efield.(pts)) / norm(Efield.(pts))
    err_ENPDL_field = norm(field_ENPDL + Efield.(pts)) / norm(Efield.(pts))

    # errors of interior problems
    @test err_IDPSL_pot < 0.01
    @test err_IDPDL_pot < 0.01
    @test err_INPSL_pot < 0.02
    @test err_INPDL_pot < 0.01

    @test err_IDPSL_field < 0.01
    @test err_IDPDL_field < 0.02
    @test err_INPSL_field < 0.02
    @test err_INPDL_field < 0.01

    # errors of exterior problems
    @test err_EDPSL_pot < 0.01
    @test err_EDPDL_pot < 0.02
    @test err_ENPSL_pot < 0.01
    @test err_ENPDL_pot < 0.01

    @test err_EDPSL_field < 0.01
    @test err_EDPDL_field < 0.03
    @test err_ENPSL_field < 0.01
    @test err_ENPDL_field < 0.02
end

## Test here some of the mixed discretizatinos
#@testset "Helmholtz potential operators: mixed discretizations with duallagrangecxd0" begin
    # r = 10.0
    # sphere = meshsphere(r, 0.2 * r)
    # λ = 20 * r
    # k = 2 * π / λ
    k = 0.031415926535897934
    sphere = readmesh(joinpath(dirname(pathof(BEAST)),"../test/assets/sphere_rad=10_h=2.in"))
    X0 = lagrangecxd0(sphere)
    X1 = lagrangec0d1(sphere)
    Y0 = duallagrangecxd0(sphere; interpolatory=true)

    S = Helmholtz3D.singlelayer(; gamma=im * k)
    D = Helmholtz3D.doublelayer(; gamma=im * k)
    Dt = Helmholtz3D.doublelayer_transposed(; gamma=im * k)
    N = Helmholtz3D.hypersingular(; gamma=im * k)

    q = 100.0
    ϵ = 1.0

    # Interior problem
    # Formulations from Sauter and Schwab, Boundary Element Methods(2011), Chapter 3.4.1.1
    r = 10.0
    pos1 = SVector(r * 1.5, 0.0, 0.0)  # positioning of point charges
    pos2 = SVector(-r * 1.5, 0.0, 0.0)

    charge1 = Helmholtz3D.monopole(position=pos1, amplitude=q/(4*π*ϵ), wavenumber=k)
    charge2 = Helmholtz3D.monopole(position=pos2, amplitude=-q/(4*π*ϵ), wavenumber=k)

    # Potential of point charges

    Φ_inc(x) = charge1(x) + charge2(x)

    gD1 = assemble(DirichletTrace(charge1), Y0) + assemble(DirichletTrace(charge2), Y0)
    gN = assemble(∂n(charge1), X1) + assemble(BEAST.n ⋅ grad(charge2), X1)

    G = assemble(Identity(), X1, Y0)
    o = ones(numfunctions(X1))

    # Interior Dirichlet problem - compare Sauter & Schwab eqs. 3.81
    M_IDPDL = (-1 / 2 * assemble(Identity(), Y0, X1) + assemble(D, Y0, X1)) # Double layer (DL)

    # Interior Neumann problem
    # Neumann derivative from SL potential with deflected nullspace
    M_INPSL = (1 / 2 * assemble(Identity(), X1, Y0) + assemble(Dt, X1, Y0)) + G * o * o' * G 

    ρ_IDPDL = M_IDPDL \ (-gD1)
    ρ_INPSL = M_INPSL \ (-gN)

    pts = meshsphere(0.8 * r, 0.8 * 0.6 * r).vertices # sphere inside on which the potential and field are evaluated
    @show length(pts)

    pot_IDPDL = potential(HH3DDoubleLayerNear(im * k), pts, ρ_IDPDL, X1; type=ComplexF64)

    pot_INPSL = potential(HH3DSingleLayerNear(im * k), pts, ρ_INPSL, Y0; type=ComplexF64)

    # Total field inside should be zero
    err_IDPDL_pot = norm(pot_IDPDL + Φ_inc.(pts)) / norm(Φ_inc.(pts))
    err_INPSL_pot = norm(pot_INPSL + Φ_inc.(pts)) / norm(Φ_inc.(pts))

    # Efield(x) = - grad Φ_inc(x)
    Efield(x) =  -grad(charge1)(x) + -grad(charge2)(x)

    field_IDPDL = -potential(HH3DHyperSingularNear(im * k), pts, ρ_IDPDL, X1)
    field_INPSL = -potential(HH3DDoubleLayerTransposedNear(im * k), pts, ρ_INPSL, Y0)

    err_IDPDL_field = norm(field_IDPDL + Efield.(pts)) / norm(Efield.(pts))
    err_INPSL_field = norm(field_INPSL + Efield.(pts)) / norm(Efield.(pts))

    # errors of interior problems
    @test err_IDPDL_pot < 0.005
    @test err_INPSL_pot < 0.002

    @test err_IDPDL_field < 0.0095
    @test err_INPSL_field < 0.002
#end

## Test here some of the mixed discretizatinos
#@testset "Helmholtz potential operators: mixed discretizations with duallagrangec0d1" begin
    # r = 10.0
    # λ = 20 * r
    # k = 2 * π / λ
    # sphere = meshsphere(r, 0.2 * r)
    k = 0.031415926535897934
    sphere = readmesh(joinpath(dirname(pathof(BEAST)),"../test/assets/sphere_rad=10_h=2.in"))

    X0 = lagrangecxd0(sphere)
    Y1 = duallagrangec0d1(sphere)

    S = Helmholtz3D.singlelayer(; gamma=im * k)
    D = Helmholtz3D.doublelayer(; gamma=im * k)
    Dt = Helmholtz3D.doublelayer_transposed(; gamma=im * k)
    N = Helmholtz3D.hypersingular(; gamma=im * k)

    q = 100.0
    ϵ = 1.0

    # Interior problem
    # Formulations from Sauter and Schwab, Boundary Element Methods(2011), Chapter 3.4.1.1
    r = 10.0
    pos1 = SVector(r * 1.5, 0.0, 0.0)  # positioning of point charges
    pos2 = SVector(-r * 1.5, 0.0, 0.0)

    charge1 = Helmholtz3D.monopole(position=pos1, amplitude=q/(4*π*ϵ), wavenumber=k)
    charge2 = Helmholtz3D.monopole(position=pos2, amplitude=-q/(4*π*ϵ), wavenumber=k)

    # Potential of point charges

    Φ_inc(x) = charge1(x) + charge2(x)

    gD1 = assemble(DirichletTrace(charge1), X0) + assemble(DirichletTrace(charge2), X0)
    gN = assemble(∂n(charge1), Y1) + assemble(BEAST.n ⋅ grad(charge2), Y1)

    G = assemble(Identity(), Y1, X0)
    o = ones(numfunctions(X0))

    # Interior Dirichlet problem - compare Sauter & Schwab eqs. 3.81
    M_IDPDL = (-1 / 2 * assemble(Identity(), X0, Y1) + assemble(D, X0, Y1)) # Double layer (DL)

    # Interior Neumann problem
    # Neumann derivative from SL potential with deflected nullspace
    M_INPSL = (1 / 2 * assemble(Identity(), Y1, X0) + assemble(Dt, Y1, X0)) + G * o * o' * G 

    ρ_IDPDL = M_IDPDL \ (-gD1)
    ρ_INPSL = M_INPSL \ (-gN)

    pts = meshsphere(0.8 * r, 0.8 * 0.6 * r).vertices # sphere inside on which the potential and field are evaluated

    pot_IDPDL = potential(HH3DDoubleLayerNear(im * k), pts, ρ_IDPDL, Y1; type=ComplexF64)

    pot_INPSL = potential(HH3DSingleLayerNear(im * k), pts, ρ_INPSL, X0; type=ComplexF64)

    # Total field inside should be zero
    err_IDPDL_pot = norm(pot_IDPDL + Φ_inc.(pts)) / norm(Φ_inc.(pts))
    err_INPSL_pot = norm(pot_INPSL + Φ_inc.(pts)) / norm(Φ_inc.(pts))

    # Efield(x) = - grad Φ_inc(x)
    Efield(x) =  -grad(charge1)(x) + -grad(charge2)(x)

    field_IDPDL = -potential(HH3DHyperSingularNear(im * k), pts, ρ_IDPDL, Y1)
    field_INPSL = -potential(HH3DDoubleLayerTransposedNear(im * k), pts, ρ_INPSL, X0)

    err_IDPDL_field = norm(field_IDPDL + Efield.(pts)) / norm(Efield.(pts))
    err_INPSL_field = norm(field_INPSL + Efield.(pts)) / norm(Efield.(pts))

    # errors of interior problems
    @test err_IDPDL_pot < 0.02
    @test err_INPSL_pot < 0.025

    @test err_IDPDL_field < 0.02
    @test err_INPSL_field < 0.025
#end