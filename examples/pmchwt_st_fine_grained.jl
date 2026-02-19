using BEAST
using BEAST.BlockArrays
using BEAST.NestedUnitRanges

using CompScienceMeshes
using LinearAlgebra
using Makeitso
using PlotlyJS
using SparseArrays

LinearAlgebra.BLAS.set_num_threads(8)


function zerofield(x)
    point(ComplexF64, 0.0im, 0.0im, 0.0im)
end
BEAST.scalartype(::typeof(zerofield)) = ComplexF64


@target geo (;h)->begin

    m1 = meshgeo(joinpath(@__DIR__, "assets/twoboxes2.geo"), h = h, physical="Gamma1")
    m2 = meshgeo(joinpath(@__DIR__, "assets/twoboxes2.geo"), h = h, physical="Gamma2")
    m3 = meshgeo(joinpath(@__DIR__, "assets/twoboxes2.geo"), h = h, physical="Gamma3")

    CompScienceMeshes.orient(m1)
    CompScienceMeshes.orient(m2)
    CompScienceMeshes.orient(m3)

    m2.vertices = m1.vertices
    m3.vertices = m1.vertices

    CompScienceMeshes.isoriented(m1)
    CompScienceMeshes.isoriented(m2)
    CompScienceMeshes.isoriented(m3)

    Σ = CompScienceMeshes.union(m1, m2, m3)

    b0 = CompScienceMeshes.union(-m1, m2)
    b1 = CompScienceMeshes.union(m1, -m3)
    b2 = CompScienceMeshes.union(-m2, m3)

    @assert CompScienceMeshes.isoriented(b0)
    @assert CompScienceMeshes.isoriented(b1)
    @assert CompScienceMeshes.isoriented(b2)

    return (;∂Ω=[b0,b1,b2], Σ)
end

@target numdoms () -> begin
    return 3
end



@target excitation (;materials, ω) -> begin

    if materials.idx != 1
        return (;Einc = zerofield, Hinc = zerofield)
    end

    κ = materials.n  * ω
    η = materials.η
    Einc = Maxwell3D.planewave(direction=(x̂+ẑ)/√2, polarization=ŷ, wavenumber=κ)
    Hinc = -1/(im*κ*η)*curl(Einc)
    return (;Einc, Hinc)
end


@target bilforms (geo, numdoms; ω, materials) -> begin

    (;idx, n, η) = materials; i = idx
    κ = n  * ω

    T = Maxwell3D.singlelayer(wavenumber=κ)
    K = Maxwell3D.doublelayer(wavenumber=κ)

    @hilbertspace m j
    @hilbertspace p q
    @hilbertspace z
    u = BEAST.hilbertspace(:u, numdoms)
    v = BEAST.hilbertspace(:v, numdoms)

    A = (
        (-K)[p,m] + η * T[p,j] -
        (1/η) * T[q,m] + (-K)[q,j]
    )[u[i],v[i]]

    (;∂Ω, Σ) = geo
    
    Edges = skeleton(Σ,1)
    edges = skeleton(∂Ω[i], 1)
    R = CompScienceMeshes.embedding(edges, Edges)
    R = R[m,p][u[i],z] + R[j,q][u[i],z]

    # @show typeof(A)

    return (;A, R)
end

@target linforms (numdoms, excitation; ω, materials) -> begin

    (;idx, n, η) = materials; i = idx
    κ = n  * ω
    (;Einc, Hinc) = excitation

    @hilbertspace p q
    v = BEAST.hilbertspace(:v, numdoms)

    n = BEAST.NormalVector()

    e = (n × Einc) × n
    h = (n × Hinc) × n

    b = -(e[p] - h[q])[v[i]]

    return (;b)
end

@target spaces (geo) -> begin

    (;∂Ω, Σ) = geo

    Edges = skeleton(Σ,1)
    edges = [skeleton(∂Ωᵢ, 1) for ∂Ωᵢ in ∂Ω]

    Nd = BEAST.nedelec(Σ, Edges)
    RT = [n×BEAST.nedelec(∂Ωᵢ, e) for (∂Ωᵢ,e) in zip(∂Ω, edges)]

    U = BEAST.DirectProductSpace([rt × rt for rt in RT])
    V = BEAST.DirectProductSpace([Nd × Nd])

    return (;U, V)
end

@sweep linmaps (spaces, !bilforms; materials=[]) -> begin

    (;A, R) = bilforms
    # @show typeof(A)
    (;U,V) = spaces

    𝗔 = assemble(A, U, U; threading=:cellcoloring)
    𝗥 = assemble(R, U, V)

    return (;𝗔, 𝗥)
end


@sweep vectors (spaces, !linforms; materials=[]) -> begin

    (;b) = linforms
    (;U,) = spaces

    𝗯 = assemble(b, U)

    return (;𝗯)
end


@target matrix (linmaps) -> begin

    𝗔 = sum(lm.𝗔 for lm in eachrow(linmaps))
    𝗥 = sum(lm.𝗥 for lm in eachrow(linmaps))

    𝐴 = Matrix(𝗔)
    𝑅 = sparse(𝗥)

    𝑆 = 𝑅' * 𝐴 * 𝑅
    ax = axes(𝗥,2)
    𝑆 = BEAST.BlockArrays.BlockedArray(𝑆, (ax,ax))
    (;𝑆, 𝑅)
end


@target solution (matrix, linmaps, vectors, spaces) -> begin

    (;𝑆, 𝑅) = matrix
    𝗥 = sum(lm.𝗥 for lm in eachrow(linmaps))

    𝗯 = sum(vt.𝗯 for vt in eachrow(vectors))
    𝑐 = 𝗥' * 𝗯

    𝑆⁻¹ = BEAST.lu(𝑆)
    𝑣 = 𝑆⁻¹ * 𝑐
    𝑢 = 𝗥 * 𝑣

    (;V, U) = spaces
    v = BEAST.FEMFunction(𝑣, V)
    u = BEAST.FEMFunction(𝑢, U)

    return (;u, v)
end


@sweep excitations (!excitation; materials=[]) -> begin
    (;excitation)
end

@target nearfield (solution, excitations; z, x, materials, ω) -> begin

    function nf(um,uj,κ,η,nts)

        Xm = um.space
        Xj = uj.space

        um = um.coeffs
        uj = uj.coeffs

        K = BEAST.MWDoubleLayerField3D(wavenumber=κ)
        T = BEAST.MWSingleLayerField3D(wavenumber=κ)

        Em = potential(K, pts, um, Xm)
        Ej = potential(T, pts, uj, Xj)

        Hm = potential(T, pts, um, Xm)
        Hj = potential(K, pts, uj, Xj)

        return -Em + η * Ej, 1/η*Hm + Hj
    end

    numdoms = length(materials)
    (;u) = solution

    (;Einc, Hinc) = excitations[3,:excitation]
    
    @hilbertspace m j
    p = BEAST.hilbertspace(:p, length(materials))

    pts = [point(x,0.5,z) for z in z, x in x]
    EH = [nf(-u[p][m], -u[p][j], ω * mat.n, mat.η, pts) for (p,mat) ∈ zip(p, materials)]
 
    for i in 1:length(materials)
        idx = excitations[i,:materials].idx
        (;Einc, Hinc) = excitations[i,:excitation]
        EH[idx][1] .-= Einc.(pts)
        EH[idx][2] .-= Hinc.(pts)
    end

    return (;E = getindex.(EH,1), H = getindex.(EH,2))
end


materials = [
    (;idx=1, n=1.0, η=1.0),
    (;idx=2, n=3.0, η=1.0),
    (;idx=3, n=4.0, η=1.0)
]

x = range(-2.0,2.0,length=150)
z = range(-1.5,2.5,length=100)
nf = make(nearfield; h=0.1, ω=2.0, x, z, materials)
exc = make(excitations; h=0.1, ω=2.0, x, z, materials)
vct = make(vectors; h=0.1, ω=2.0, x, z, materials)

Etot = sum(nf.E)
Htot = sum(nf.H)

using LinearAlgebra
using PlotlyJS
hm1 = PlotlyJS.heatmap(x=x, y=z, z=real.(getindex.(nf.E[1],2)), colorscale="Viridis", zmin=-2, zmax=2)
hm2 = PlotlyJS.heatmap(x=x, y=z, z=real.(getindex.(nf.E[2],2)), colorscale="Viridis", zmin=-2, zmax=2)
hm3 = PlotlyJS.heatmap(x=x, y=z, z=real.(getindex.(nf.E[3],2)), colorscale="Viridis", zmin=-2, zmax=2)
hm4 = PlotlyJS.heatmap(x=x, y=z, z=real.(getindex.(Etot,2)),  colorscale="Viridis", zmin=-2, zmax=2)

plt = Plot(Layout(Subplots(rows=2,cols=2, specs=[Spec() Spec(); Spec() Spec()])));
PlotlyJS.add_trace!(plt, hm1, row=1, col=1);
PlotlyJS.add_trace!(plt, hm2, row=1, col=2);
PlotlyJS.add_trace!(plt, hm3, row=2, col=1);
PlotlyJS.add_trace!(plt, hm4, row=2, col=2);
display(plt)
