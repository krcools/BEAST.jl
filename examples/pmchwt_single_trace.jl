using BEAST
using BEAST.BlockArrays
using BEAST.NestedUnitRanges

using CompScienceMeshes
using LinearAlgebra
using Makeitso
using PlotlyJS
using SparseArrays

LinearAlgebra.BLAS.set_num_threads(8)

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


@target material (;κ) -> begin
    return [
        (;κ=κ, η=1.0),
        (;κ=2*κ, η=1.0),
        (;κ=3*κ, η=1.0)
    ]
end

@target excitation (material) -> begin
    Einc = Maxwell3D.planewave(direction=(x̂+ẑ)/√2, polarization=ŷ, wavenumber=material[1].κ)
    Hinc = -1/(im*material[1].κ*material[1].η)*curl(Einc)
    return [
        (;Einc=Einc, Hinc=Hinc),
        (;Einc=nothing, Hinc=nothing),
        (;Einc=nothing, Hinc=nothing)
    ]
end


@target formulation (geo, material, excitation) -> begin

    numdoms = length(material)
    (;∂Ω, Σ) = geo

    (;κ, η) = first(material)
    (;Einc, Hinc) = first(excitation)
    
    e = (n × Einc) × n
    h = (n × Hinc) × n

    Ts = [Maxwell3D.singlelayer(wavenumber=m.κ) for m in material]
    Ks = [Maxwell3D.doublelayer(wavenumber=m.κ) for m in material]
    
    @hilbertspace m j
    @hilbertspace p q
    @hilbertspace z
    u = BEAST.hilbertspace(:u, numdoms)
    v = BEAST.hilbertspace(:v, numdoms)

    b = -(e[p] - h[q])[v[1]]

    Aloc = [
        (-K)[p,m] + mat.η * T[p,j] -
        (1/mat.η) * T[q,m] + (-K)[q,j] for (T,K,mat) in zip(Ts,Ks,material)
    ]
    A = sum(Aloc[i][u[i],v[i]] for i in 1:numdoms)

    Edges = skeleton(Σ,1)
    edges = [skeleton(∂Ωᵢ, 1) for ∂Ωᵢ in ∂Ω]
    Rloc = [ CompScienceMeshes.embedding(edges[i], Edges) for i  in 1:numdoms]
    R = sum((Rloc[i][m,p])[u[i],z] + (Rloc[i][j,q])[u[i],z] for i in 1:numdoms) 

    return (;bilforms=(;A,R), linforms=(;b))
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


@target discretization (formulation, spaces) -> begin

    (;bilforms, linforms) = formulation
    (;A,R) = bilforms
    (;b) = linforms

    (;U,V) = spaces

    bu = assemble(b, U)
    RUv = assemble(R, U, V)
    Auu = assemble(A, U, U)

    matrices = (;Auu, RUv)
    vectors = (;bu)
    return (;matrices, vectors)
end


@target solution (discretization, spaces) -> begin

    (;matrices, vectors) = discretization
    (;Auu, RUv) = matrices
    (;bu) = vectors

    (;U,V) = spaces

    MA = Matrix(Auu)
    MR = sparse(RUv)

    A = MR' * MA * MR
    A = BEAST.BlockArrays.BlockedArray(A, (
        BEAST.NestedUnitRanges.nestedrange(V, 1, numfunctions) ,
        BEAST.NestedUnitRanges.nestedrange(V, 1, numfunctions) ) )
    b = RUv' * bu

    Ai = BEAST.lu(A)
    v = Ai * b
    u = RUv * v

    v = BEAST.FEMFunction(v, V)
    u = BEAST.FEMFunction(u, U)

    return (;u)
end

@target nearfield (solution, material, excitation; z, x) -> begin

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

    numdoms = length(material)
    (;u) = solution
    (;Einc, Hinc) = first(excitation)
    
    @hilbertspace m j
    p = BEAST.hilbertspace(:p, length(material))

    pts = [point(x,0.5,z) for z in z, x in x]
    EH = [nf(-u[p][m], -u[p][j], mat.κ, mat.η, pts) for (p,mat) ∈ zip(p, material)]
 
    EH[1][1] .-= Einc.(pts)
    EH[1][2] .-= Hinc.(pts)

    return (;E = getindex.(EH,1), H = getindex.(EH,2))
end


x = range(-2.0,2.0,length=150)
z = range(-1.5,2.5,length=100)
nf = make(nearfield; h=0.08, κ=2.0, x=x, z=z)

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


## Debugging

# ds = make(discretization; h=0.08, κ=2.0)
ds = make(discretization; h=0.15, κ=2.0)

using LinearAlgebra
using BEAST.LinearMaps
using BEAST.BlockArrays
using BEAST.BlockArrays.ArrayLayouts
using SparseArrays

(;matrices, vectors) = ds
(;Auu, RUv) = matrices
(;bu) = vectors

MA1 = AbstractMatrix(Auu)
MR1 = AbstractMatrix(RUv)
# vb = AbstractArray(bu)

MA2 = Matrix(Auu)
MR2 = sparse(RUv)

MulAB = ArrayLayouts.Mul(MA1, MR1)
@which ArrayLayouts.materialize(MulAB)
InstAB = ArrayLayouts.instantiate(MulAB)
@which copy(InstAB)
@which ArrayLayouts.mulreduce(InstAB)
redAB = ArrayLayouts.mulreduce(InstAB)
@which copy(redAB)
simAB = similar(redAB)
@which copyto!(simAB, redAB)

M = redAB
dest = simAB
fct = ArrayLayouts._fill_copyto!(dest, M.C)
@which muladd!(M.α, M.A, M.B, M.β, fct; Czero = M.Czero)
error()

@enter copyto!(simAB, redAB)

@time MA2 * MR2
@time MA1 * MR1

axs1 = axes(Auu)
axs2 = Tuple( blockedrange([sz[1] for sz in blocksizes(ax1)]) for ax1 in axs1 )

T = ComplexF64
Y1 = similar(Array{T}, axs1)
Y2 = similar(Array{T}, axs2)

@which LinearAlgebra.rmul!(Y1, true)
@which LinearAlgebra.rmul!(Y2, true)

@profview LinearAlgebra.rmul!(Y1, true)
@profview LinearAlgebra.rmul!(Y2, true)

lmap = Auu.maps[1].lmap
P = getindex(axes(lmap,1), lmap.I)
Q = getindex(axes(lmap,2), lmap.J)

Y1IJ = view(Y1, P, Q)
Y2IJ = view(Y2, P, Q)

@which LinearAlgebra.rmul!(Y1IJ, true)
@which LinearAlgebra.rmul!(Y2IJ, true)

@profview LinearAlgebra.rmul!(Y1IJ, true)
@profview LinearAlgebra.rmul!(Y2IJ, true)

# error()
@time LinearMaps._unsafe_mul!(Y1, Auu, false)
@profview LinearMaps._unsafe_mul!(Y1, Auu, false)
# @profview AbstractMatrix{T}(Auu)

error()