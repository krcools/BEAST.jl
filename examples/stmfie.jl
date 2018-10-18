using CompScienceMeshes, BEAST
o, x, y, z = euclidianbasis(3)

# D, Δx = 1.0, 0.35

Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
X, Y = raviartthomas(Γ), buffachristiansen(Γ)

# Δt, Nt = 0.11, 200
Δt, Nt = 0.6, 200
δ = timebasisdelta(Δt, Nt)
h = timebasisc0d1(Δt, Nt)
p = timebasiscxd0(Δt, Nt)
V = X ⊗ p
W = Y ⊗ p

width, delay, scaling = 8.0, 12.0, 1.0
gaussian = creategaussian(width, delay, scaling)

direction, polarisation = ẑ, x̂
E = BEAST.planewave(polarisation, direction, gaussian, 1.0)
# E = BEAST.planewave(polarisation, direction, derive(gaussian), 1.0)
H = direction × E

@hilbertspace j; @hilbertspace m′
K, I, N = MWDoubleLayerTDIO(1.0, 1.0, 0), Identity(), NCross()
tdmfie = @discretise (0.5*(N⊗I) + 1.0*K)[m′,j] == -1H[m′]   j∈V  m′∈W
xmfie = solve(tdmfie)

Xmfie, Δω, ω0 = fouriertransform(xmfie, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω .- 1.0))

ω1 = ω[i1]
um = Xmfie[:,i1] / fouriertransform(gaussian)(ω1)
