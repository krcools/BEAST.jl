using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Plotly

## Looking at convergence
hs = [0.3, 0.2, 0.1, 0.09]#,0.08,0.07,0.06,0.04]
ir = 0.8
err_IDPSL_pot = zeros(Float64, length(hs))
err_IDPDL_pot = zeros(Float64, length(hs))
err_INPSL_pot = zeros(Float64, length(hs))
err_INPDL_pot = zeros(Float64, length(hs))
err_IDPSL_field = zeros(Float64, length(hs))
err_IDPDL_field = zeros(Float64, length(hs))
err_INPSL_field = zeros(Float64, length(hs))
err_INPDL_field = zeros(Float64, length(hs))

err_EDPSL_pot = zeros(Float64, length(hs))
err_EDPDL_pot = zeros(Float64, length(hs))
err_ENPSL_pot = zeros(Float64, length(hs))
err_ENPDL_pot = zeros(Float64, length(hs))

err_EDPSL_field = zeros(Float64, length(hs))
err_EDPDL_field = zeros(Float64, length(hs))
err_ENPSL_field = zeros(Float64, length(hs))
err_ENPDL_field = zeros(Float64, length(hs))

for (i, h) in enumerate(hs)
    r = 50.0
    sphere = meshsphere(r, h * r)
    X0 = lagrangecxd0(sphere)
    X1 = lagrangec0d1(sphere)

    S = Helmholtz3D.singlelayer(; gamma=0.0)
    D = Helmholtz3D.doublelayer(; gamma=0.0)
    Dt = Helmholtz3D.doublelayer_transposed(; gamma=0.0)
    N = Helmholtz3D.hypersingular(; gamma=0.0)

    q = 100.0
    ϵ = 1.0
# Interior problem
# Formulations from Sauter and Schwab, Boundary Element Methods(2011), Chapter 3.4.1.1

    pos1 = SVector(r * 1.5, 0.0, 0.0)
    pos2 = SVector(-r * 1.5, 0.0, 0.0)

    charge1 = HH3DMonopole(position=pos1, amplitude=q/(4*π*ϵ), wavenumber=0.0)
    charge2 = HH3DMonopole(position=pos2, amplitude=-q/(4*π*ϵ), wavenumber=0.0)

    #Potential of point charges
    Φ_inc(x) = charge1(x) + charge2(x)

    gD0 = assemble(DirichletTrace(charge1), X0) + assemble(DirichletTrace(charge2), X0)
    gD1 = assemble(DirichletTrace(charge1), X1) + assemble(DirichletTrace(charge2), X1)
    gN = assemble(∂n(charge1), X1) + assemble(BEAST.n ⋅ grad(charge2), X1)

    G = assemble(Identity(), X1, X1)
    o = ones(numfunctions(X1))
    #Interior Dirichlet problem - compare Sauter & Schwab
    M_IDPSL = assemble(S, X0, X0)
    M_IDPDL = (-1 / 2 * assemble(Identity(), X1, X1) + assemble(D, X1, X1))

    #Interior Neumann problem
    M_INPSL = (1 / 2 * assemble(Identity(), X1, X1) + assemble(Dt, X1, X1)) + G * o * o' * G
    M_INPDL = assemble(N, X1, X1) + G * o * o' * G

    ρ_IDPSL = M_IDPSL \ (-gD0)
    ρ_IDPDL = M_IDPDL \ (-gD1)
    ρ_INPSL = M_INPSL \ (-gN)
    ρ_INPDL = M_INPDL \ (gN)

    pts = meshsphere(r * ir, r * ir * 0.6).vertices

    pot_IDPSL = potential(HH3DSingleLayerNear(0.0), pts, ρ_IDPSL, X0; type=ComplexF64)
    pot_IDPDL = potential(HH3DDoubleLayerNear(0.0), pts, ρ_IDPDL, X1; type=ComplexF64)
    pot_INPSL = potential(HH3DSingleLayerNear(0.0), pts, ρ_INPSL, X1; type=ComplexF64)
    pot_INPDL = potential(HH3DDoubleLayerNear(0.0), pts, ρ_INPDL, X1; type=ComplexF64)

    # Total potential inside should be zero
    err_IDPSL_pot[i] = norm(pot_IDPSL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))
    err_IDPDL_pot[i] = norm(pot_IDPDL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))
    err_INPSL_pot[i] = norm(pot_INPSL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))
    err_INPDL_pot[i] = norm(pot_INPDL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))

    Efield(x) =  -grad(charge1)(x) + -grad(charge2)(x)

    field_IDPSL = -potential(HH3DDoubleLayerTransposedNear(0.0), pts, ρ_IDPSL, X0)
    field_IDPDL = -potential(HH3DHyperSingularNear(0.0), pts, ρ_IDPDL, X1)
    field_INPSL = -potential(HH3DDoubleLayerTransposedNear(0.0), pts, ρ_INPSL, X1)
    field_INPDL = -potential(HH3DHyperSingularNear(0.0), pts, ρ_INPDL, X1)

    # Total field inside should be zero
    err_IDPSL_field[i] = norm(field_IDPSL + Efield.(pts)) / norm(Efield.(pts))
    err_IDPDL_field[i] = norm(field_IDPDL + Efield.(pts)) / norm(Efield.(pts))
    err_INPSL_field[i] = norm(field_INPSL + Efield.(pts)) / norm(Efield.(pts))
    err_INPDL_field[i] = norm(field_INPDL + Efield.(pts)) / norm(Efield.(pts))

    # Exterior problem
    # formulations from Sauter and Schwab, Boundary Element Methods(2011), Chapter 3.4.1.2

    pos1 = SVector(r * 0.5, 0.0, 0.0)
    pos2 = SVector(-r * 0.5, 0.0, 0.0)

    charge1 = HH3DMonopole(position=pos1, amplitude=q/(4*π*ϵ), wavenumber=0.0)
    charge2 = HH3DMonopole(position=pos2, amplitude=-q/(4*π*ϵ), wavenumber=0.0)

    gD0 = assemble(DirichletTrace(charge1), X0) + assemble(DirichletTrace(charge2), X0)
    gD1 = assemble(DirichletTrace(charge1), X1) + assemble(DirichletTrace(charge2), X1)
    gN = assemble(∂n(charge1), X1) + assemble(∂n(charge2), X1)

    G = assemble(Identity(), X1, X1)
    o = ones(numfunctions(X1))

    M_EDPSL = assemble(S, X0, X0)
    M_EDPDL = (1 / 2 * assemble(Identity(), X1, X1) + assemble(D, X1, X1))

    M_ENPSL = (-1 / 2 * assemble(Identity(), X1, X1) + assemble(Dt, X1, X1)) + G * o * o' * G
    M_ENPDL = assemble(N, X1, X1) + G * o * o' * G

    ρ_EDPSL = M_EDPSL \ (-gD0)
    ρ_EDPDL = M_EDPDL \ (-gD1)
    ρ_ENPSL = M_ENPSL \ (-gN)
    ρ_ENPDL = M_ENPDL \ (gN)

    testsphere = meshsphere(r / ir, r / ir * 0.6)
    pts = testsphere.vertices[norm.(testsphere.vertices) .> r]

    pot_EDPSL = potential(HH3DSingleLayerNear(0.0), pts, ρ_EDPSL, X0; type=ComplexF64)
    pot_EDPDL = potential(HH3DDoubleLayerNear(0.0), pts, ρ_EDPDL, X1; type=ComplexF64)
    pot_ENPSL = potential(HH3DSingleLayerNear(0.0), pts, ρ_ENPSL, X1; type=ComplexF64)
    pot_ENPDL = potential(HH3DDoubleLayerNear(0.0), pts, ρ_ENPDL, X1; type=ComplexF64)

    err_EDPSL_pot[i] = norm(pot_EDPSL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))
    err_EDPDL_pot[i] = norm(pot_EDPDL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))
    err_ENPSL_pot[i] = norm(pot_ENPSL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))
    err_ENPDL_pot[i] = norm(pot_ENPDL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))

    field_EDPSL = -potential(HH3DDoubleLayerTransposedNear(0.0), pts, ρ_EDPSL, X0)
    field_EDPDL = -potential(HH3DHyperSingularNear(0.0), pts, ρ_EDPDL, X1)
    field_ENPSL = -potential(HH3DDoubleLayerTransposedNear(0.0), pts, ρ_ENPSL, X1)
    field_ENPDL = -potential(HH3DHyperSingularNear(0.0), pts, ρ_ENPDL, X1)

    err_EDPSL_field[i] = norm(field_EDPSL + Efield.(pts)) / norm(Efield.(pts))
    err_EDPDL_field[i] = norm(field_EDPDL + Efield.(pts)) / norm(Efield.(pts))
    err_ENPSL_field[i] = norm(field_ENPSL + Efield.(pts)) / norm(Efield.(pts))
    err_ENPDL_field[i] = norm(field_ENPDL + Efield.(pts)) / norm(Efield.(pts))
end

##
Plotly.plot(
    [
        Plotly.scatter(; x=hs, y=err_IDPSL_pot[:, end], name="IDPSL"),
        Plotly.scatter(; x=hs, y=err_IDPDL_pot[:, end], name="IDPDL"),
        Plotly.scatter(; x=hs, y=err_INPSL_pot[:, end], name="INPSL"),
        Plotly.scatter(; x=hs, y=err_INPDL_pot[:, end], name="INPDL"),
    ],
    Layout(;
        xaxis_type="log",
        yaxis_type="log",
        xaxis_title="h / r",
        yaxis_title="rel. error",
        title="Errors - potential - interior problem",
    ),
)

Plotly.plot(
    [
        Plotly.scatter(; x=hs, y=err_EDPSL_pot[:, end], name="EDPSL"),
        Plotly.scatter(; x=hs, y=err_EDPDL_pot[:, end], name="EDPDL"),
        Plotly.scatter(; x=hs, y=err_ENPSL_pot[:, end], name="ENPSL"),
        Plotly.scatter(; x=hs, y=err_ENPDL_pot[:, end], name="ENPDL"),
    ],
    Layout(;
        xaxis_type="log",
        yaxis_type="log",
        xaxis_title="h / r",
        yaxis_title="rel. error",
        title="Errors - potential - exterior problem",
    ),
)

Plotly.plot(
    [
        Plotly.scatter(; x=hs, y=err_IDPSL_field[:, end], name="IDPSL"),
        Plotly.scatter(; x=hs, y=err_IDPDL_field[:, end], name="IDPDL"),
        Plotly.scatter(; x=hs, y=err_INPSL_field[:, end], name="INPSL"),
        Plotly.scatter(; x=hs, y=err_INPDL_field[:, end], name="INPDL"),
    ],
    Layout(;
        xaxis_type="log",
        yaxis_type="log",
        xaxis_title="h / r",
        yaxis_title="rel. error",
        title="Errors - field - interior problem",
    ),
)

Plotly.plot(
    [
        Plotly.scatter(; x=hs, y=err_EDPSL_field[:, end], name="EDPSL"),
        Plotly.scatter(; x=hs, y=err_EDPDL_field[:, end], name="EDPDL"),
        Plotly.scatter(; x=hs, y=err_ENPSL_field[:, end], name="ENPSL"),
        Plotly.scatter(; x=hs, y=err_ENPDL_field[:, end], name="ENPDL"),
    ],
    Layout(;
        xaxis_type="log",
        yaxis_type="log",
        xaxis_title="h / r",
        yaxis_title="rel. error",
        title="Errors - field - exterior problem",
    ),
)