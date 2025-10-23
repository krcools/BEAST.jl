using BEAST
using CompScienceMeshes
using Makeitso
using PlotlyJS

@target geo (;h) -> begin
    (;∂Ω = [meshcuboid(1.0, 1.0, 1.0, h)])
end

@target material (;κ) -> begin
    return [
        (;κ=κ, η=1.0),
        (;κ=3κ, η=1.0),
    ]
end

@target excitation (material) -> begin

    κ = material[1].κ
    Uinc = Helmholtz3D.planewave(wavenumber=κ, direction=point(0,0,1))
    return [(;Uinc)]
end

@target formulation (material, excitation, geo) -> begin

    numdoms = 2
    Γ = geo.∂Ω[1]

    # 𝐼 = BEAST.Identity()
    𝑆 = [Helmholtz3D.singlelayer(wavenumber=m.κ) for m in material]
    𝐾 = [Helmholtz3D.doublelayer_transposed(wavenumber=m.κ) for m in material]
    𝐷 = [Helmholtz3D.doublelayer(wavenumber=m.κ) for m in material]
    𝑊 = [Helmholtz3D.hypersingular(wavenumber=m.κ) for m in material]

    @hilbertspace u v
    @hilbertspace p q

    uinc = strace(excitation[1].Uinc, Γ)
    vinc = ∂n(excitation[1].Uinc)

    𝐴 = [
        𝐷[i][p,u] + 𝑆[i][p,v] +
        𝑊[i][q,u] - 𝐾[i][q,v] for i in 1:numdoms]

    𝑀 = 𝐴[1] + 𝐴[2]
    𝑏 = -uinc[p] + vinc[q]

    return (;bilforms=(;𝑀), linforms=(;𝑏))
end

@target spaces (geo) -> begin
    (;∂Ω) = geo
    Γ = ∂Ω[1]

    L = lagrangec0d1(Γ)
    P = lagrangecxd0(Γ)

    U = L × P
    V = P × L
    return (;U=U, V=V)
end

@target discretization (formulation, spaces) -> begin
    (;bilforms, linforms) = formulation
    (;𝑀) = bilforms
    (;𝑏) = linforms

    (;U, V) = spaces

    M = assemble(𝑀, V, U)
    b = assemble(𝑏, V)

    return (;matrices=(;M), vectors=(;b))
end

@target solution (discretization, spaces) -> begin
    (;matrices, vectors) = discretization
    (;M) = matrices
    (;b) = vectors

    (;U, V) = spaces

    𝗠 = Matrix(M)
    𝗯 = Vector(b)

    𝘂 = 𝗠 \ 𝗯
    u = BEAST.BlockArrays.BlockedVector(𝘂, (
        BEAST.NestedUnitRanges.nestedrange(U, 1, numfunctions),))
    u = BEAST.FEMFunction(u, U)
    return (;u)
end

@target nearfield (solution, material, excitation; pts) -> begin

    function nf(um,uj,mat,pts)

        (;κ) = mat

        Xm = um.space
        Xj = uj.space

        um = um.coeffs
        uj = uj.coeffs

        𝒟 = BEAST.HH3DDoubleLayerNear(wavenumber=κ)
        𝒮 = BEAST.HH3DSingleLayerNear(wavenumber=κ)

        U1 = potential(𝒟, pts, um, Xm; type=ComplexF64)
        U2 = potential(𝒮, pts, uj, Xj; type=ComplexF64)

        return U1 + U2
    end

    (;u) = solution
    (;Uinc) = excitation[1]

    @hilbertspace p v
    U = [
        nf(u[p], u[v], material[1], pts),
        nf(-u[p], -u[v], material[2], pts),]
    U[1] .+= Uinc.(pts)

    return (;U)
end

x = range(-2.0,2.0,length=75)
z = range(-1.5,2.5,length=50)
pts = [point(xi, 0.5, zi) for xi in x, zi in z]
nf = make(nearfield; h=0.1, κ=3.0, pts=pts)

using LinearAlgebra
using PlotlyJS
hm1 = PlotlyJS.heatmap(x=x, y=z, z=real.(nf.U[1]), colorscale="Viridis", zmin=-2, zmax=2)
hm2 = PlotlyJS.heatmap(x=x, y=z, z=real.(nf.U[2]), colorscale="Viridis", zmin=-2, zmax=2)

plt = Plot(Layout(Subplots(rows=2,cols=1, specs=reshape([Spec(); Spec()], 2, 1))));
PlotlyJS.add_trace!(plt, hm1, row=1, col=1);
PlotlyJS.add_trace!(plt, hm2, row=2, col=1);
display(plt)