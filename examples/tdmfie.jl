using PlotlyJS
include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))

using CompScienceMeshes, BEAST
o, x, y, z = euclidianbasis(3)

#D, Δx = 1.0, 0.25
D, Δx = 1.0, 0.35
Γ = meshsphere(D, Δx)
X, Y = raviartthomas(Γ), buffachristiansen(Γ)

#Δt, Nt = 0.1, 600
Δt, Nt = 0.11, 200
T, δ = timebasisshiftedlagrange(Δt, Nt, 2), timebasisdelta(Δt, Nt)

V = X ⊗ T; W = Y ⊗ δ

width, delay, scaling = 8.0, 12.0, 1.0
gaussian = creategaussian(width, delay, scaling)

direction, polarisation = z, x
E = planewave(polarisation, direction, gaussian, 1.0); H = direction × E

@hilbertspace j; @hilbertspace m′
K, I, N = MWDoubleLayerTDIO(1.0, 1.0, 0), Identity(), NCross()
tdmfie = @discretise (0.5*(N⊗I) + 1.0*K)[m′,j] == -1E[m′]   j∈V  m′∈W
xmfie = solve(tdmfie)


t3 = scatter(x=0:Δt:(Nt-1)Δt, y=xmfie[1,:])

Xmfie, Δω, ω0 = fouriertransform(xmfie, Δt, 0.0, 2)
ω = collect(ω0 + (0:Nt-1)*Δω)
_, i1 = findmin(abs(ω-1.0))

ω1 = ω[i1]
um = Xmfie[:,i1] / fouriertransform(gaussian)(ω1)

fcrm, geo = facecurrents(um, X)
t4 = patch(geo, real.(norm.(fcrm)))
PlotlyJS.plot(t4)
