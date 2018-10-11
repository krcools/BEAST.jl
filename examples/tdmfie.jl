using CompScienceMeshes, BEAST
o, x, y, z = euclidianbasis(3)

sol = 5.0
Δt, Nt = 0.12/sol,400

D, Δx = 1.0, 0.25
Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
X, Y = raviartthomas(Γ), buffachristiansen(Γ)

T = timebasisshiftedlagrange(Δt, Nt, 2)
δ = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = Y ⊗ δ

width, delay, scaling = 8.0/sol, 12.0/sol, 1.0
gaussian = creategaussian(width, delay, scaling)

direction, polarisation = z, x
# E = planewave(polarisation, direction, gaussian, 1.0)
E = BEAST.planewave(polarisation, direction, derive(gaussian), sol)
H = direction × E

@hilbertspace j; @hilbertspace m′
K, I, N = MWDoubleLayerTDIO(sol, 1.0, 0), Identity(), NCross()
tdmfie = @discretise (0.5*(N⊗I) + 1.0*K)[m′,j] == -1H[m′]   j∈V  m′∈W
xmfie = solve(tdmfie)


using Plots
scatter(0:Δt:(Nt-1)Δt, xmfie[1,:])

Xmfie, Δω, ω0 = fouriertransform(xmfie, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω .- 1.0))

ω1 = ω[i1]
um = Xmfie[:,i1] / fouriertransform(gaussian)(ω1)
