using CompScienceMeshes
using BEAST

Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
Γref = CompScienceMeshes.barycentric_refinement(Γ)
X = raviartthomas(Γ, BEAST.Continuity{:none})
Y = raviartthomas(Γref, BEAST.Continuity{:none})

# X = subset(X,1690:1692)
# Y = subset(Y,10141:10143)

# X = subset(X,[1692])
# Y = subset(Y,[10147])

# Y = buffachristiansen(Γ)

κ, η = 1.0, 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

g = 6
# ssstrat(g) = BEAST.DoubleNumSauterQstrat(7, 6, g, g, g, g)
ssstrat(g) = BEAST.CommonFaceVertexSauterCommonEdgeWiltonPostitiveDistanceNumQStrat(
    7, 6, 10, 8, g, g, g+3, g)

qs1 = ssstrat(g)
qs2 = BEAST.NonConformingIntegralOpQStrat(ssstrat(g))

A1 = assemble(t,Y,X, quadstrat=qs1, threading=BEAST.Threading{:single})
A2 = assemble(t,Y,X, quadstrat=qs2, threading=BEAST.Threading{:single})

# @enter assemble(t,Y,X, quadstrat=qs2, threading=BEAST.Threading{:single})

import Plots
Plots.plotly()
step = 1
rowidx = 10894
Plots.plot(imag.(A1[rowidx,1:step:end]))
Plots.scatter!(imag.(A2[rowidx,1:step:end]))
val, pos = findmax(abs.(A1-A2))

rowidx = Y.fns[pos[1]][1].cellid
colidx = X.fns[pos[2]][1].cellid

τ = chart(Γref, rowidx)
σ = chart(Γ, colidx)

# τ = [
# [0.18517788982868066, 0.05998671122432393, 0.9735312677690339],
# [0.1647116144395066, 0.0050896055947335095, 0.9807011064756881],
# [0.245326672427806, -0.06121735585211891, 0.9675056894602609]]

# σ = [
# [0.1133762981558132, -0.1176791567283826, 0.9865583769286949],
# [0.1647116144395066, 0.0050896055947335025, 0.9807011064756881],
# [0.08409655645120719, 0.07139656704158594, 0.9938965234911153]]

# eτ = simplex(
#     point(0.1647116144395066, 0.0050896055947335095, 0.9807011064756881),
#     point(0.245326672427806, -0.06121735585211891, 0.9675056894602609))

# eσ = simplex(
#     point(0.1647116144395066, 0.0050896055947335025, 0.9807011064756881),
#     point(0.08409655645120719, 0.07139656704158594, 0.9938965234911153))

# CompScienceMeshes.overlap(eτ, eσ)

# τ = τ.vertices
# σ = σ.vertices

import Plotly
function Plotly.mesh3d(a::Vector{<:CompScienceMeshes.Simplex}; opacity=0.5, kwargs...)
    T = coordtype(a[1])
    n = length(a)
    X = Vector{T}(undef, 3*n)
    Y = Vector{T}(undef, 3*n)
    Z = Vector{T}(undef, 3*n)
    for i in 1:n
        X[3*(i-1)+1:3*i] = getindex.(a[i].vertices, 1)
        Y[3*(i-1)+1:3*i] = getindex.(a[i].vertices, 2)
        Z[3*(i-1)+1:3*i] = getindex.(a[i].vertices, 3)
    end

    I = collect(0:3:3*(n-1))
    J = I .+ 1
    K = I .+ 2

    return Plotly.mesh3d(x=X,y=Y,z=Z,i=I,j=J,k=K; opacity, kwargs...)
end

m1 = Plotly.mesh3d([τ], opacity=0.5)
m2 = Plotly.mesh3d([σ], opacity=0.5)


# m1 = Plotly.mesh3d(
#     x=getindex.(τ,1),
#     y=getindex.(τ,2),
#     z=getindex.(τ,3),
#     i=[0], j=[1], k=[2])
# m2 = Plotly.mesh3d(
#     x=getindex.(σ,1),
#     y=getindex.(σ,2),
#     z=getindex.(σ,3),
#     i=[0], j=[1], k=[2])
Plotly.plot([m1,m2])




isct1, clps1 = CompScienceMeshes.intersection_keep_clippings(τ,σ)
isct2, clps2 = CompScienceMeshes.intersection_keep_clippings(σ,τ)


m3 = Plotly.mesh3d(isct1, opacity=0.5, color="blue")
m4 = Plotly.mesh3d(isct2, opacity=0.5, color="red")
Plotly.plot([m3,m4])

m5 = Plotly.mesh3d(clps2[1], opacity=0.5, color="green")
m6 = Plotly.mesh3d(clps2[1][[1]], opacity=0.5, color="yellow")
Plotly.plot([m3,m6])

p = isct1[1]
q = clps2[1][1]
for (i,λ) in pairs(faces(p))
    for (j,μ) in pairs(faces(q))
        if CompScienceMeshes.overlap(λ, μ)
            global gi = i
            global gj = j
            global qr = BEAST.NonConformingTouchQRule(qs1, i, j)
end end end

𝒳 = refspace(X)
𝒴 = refspace(Y)
out10 = zeros(ComplexF64, 3, 3)
BEAST.momintegrals!(t, 𝒴, 𝒳, p, q, out10, BEAST.NonConformingTouchQRule(ssstrat(10), gi, gj))
out20 = zeros(ComplexF64, 3, 3)
BEAST.momintegrals!(t, 𝒴, 𝒳, p, q, out20, BEAST.NonConformingTouchQRule(ssstrat(20), gi, gj))

@show out10[1,1] out20[1,1]

τs, σs = BEAST._conforming_refinement_touching_triangles(p,q,2,3)
@assert length(τs) == 1
@assert length(σs) == 2
@assert volume(p) ≈ sum(volume.(τs))
@assert volume(q) ≈ sum(volume.(σs))
@assert all(volume.(τs) .> eps(Float64)*1e3)
@assert all(volume.(σs) .> eps(Float64)*1e3)
m7 = Plotly.mesh3d(τs[[1]], color="blue", opacity=0.5)
m8 = Plotly.mesh3d(σs[[1]], color="red", opacity=0.5)
Plotly.plot([m7,m8])

out1_10 = zeros(ComplexF64, 3, 3)
qs = ssstrat(10)
qd = quaddata(t, 𝒴, 𝒳, τs, σs, qs)
qr = quadrule(t, 𝒴, 𝒳, 1, τs[1], 1, σs[1], qd, qs)
BEAST.momintegrals!(t, 𝒴, 𝒳, τs[1], σs[1], out1_10, qr)
out1_20 = zeros(ComplexF64, 3, 3)
qs = ssstrat(20)
qd = quaddata(t, 𝒴, 𝒳, τs, σs, qs)
qr = quadrule(t, 𝒴, 𝒳, 1, τs[1], 1, σs[1], qd, qs)
BEAST.momintegrals!(t, 𝒴, 𝒳, τs[1], σs[1], out1_20, qr)
@show out1_10[1,1] out1_20[1,1]

out2_10 = zeros(ComplexF64, 3, 3)
qs = ssstrat(10)
qd = quaddata(t, 𝒴, 𝒳, τs, σs, qs)
qr = quadrule(t, 𝒴, 𝒳, 1, τs[1], 1, σs[2], qd, qs)
BEAST.momintegrals!(t, 𝒴, 𝒳, τs[1], σs[2], out2_10, qr)
out2_20 = zeros(ComplexF64, 3, 3)
qs = ssstrat(20)
qd = quaddata(t, 𝒴, 𝒳, τs, σs, qs)
qr = quadrule(t, 𝒴, 𝒳, 1, τs[1], 2, σs[2], qd, qs)
BEAST.momintegrals!(t, 𝒴, 𝒳, τs[1], σs[2], out2_20, qr)
@show out2_10[1,1] out2_20[1,1]


out_ref = zeros(ComplexF64, 3, 3)
qrss = BEAST.SauterSchwabQuadrature.CommonVertex(g)
BEAST.momintegrals_test_refines_trial!(out_ref, t,
    Y, rowidx, p,
    X, colidx, q,
    qrss, qs1)
@show out_ref