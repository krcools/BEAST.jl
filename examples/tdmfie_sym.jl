using CompScienceMeshes, BEAST
Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
# Γ = meshsphere(1.0, 0.4)

X = raviartthomas(Γ)
Y = BEAST.nedelec(Γ)

Δt, Nt = 0.3, 200
T = timebasisshiftedlagrange(Δt, Nt, 2)
δ = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = Y ⊗ δ

duration = 20 * Δt * 2
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)

direction, polarisation = ẑ, x̂
E = BEAST.planewave(polarisation, direction, gaussian, 1.0)
H = direction × E

DL = TDMaxwell3D.doublelayer(speedoflight=1.0)
I = Identity()
N = NCross()

@hilbertspace k
@hilbertspace j
mfie = @discretise (0.5(N⊗I) + 1.0DL)[k,j] == -1.0H[k] k∈W j∈V

# Kyx = assemble(DL, W, V)
# Nyx = assemble(N⊗I, W, V)

xmfie_sym = solve(mfie)
# xmfie_sym = xmfie

Xmfie_sym, Δω, ω0 = fouriertransform(xmfie_sym, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω .- 1.0))

ω1 = ω[i1]
um = Xmfie_sym[:,i1] / fouriertransform(gaussian)(ω1)
