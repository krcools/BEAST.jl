using CompScienceMeshes, BEAST

# Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
# Γ = meshsphere(radius=1.0, h=0.4)
X = raviartthomas(Γ)

Δt, Nt = 0.1, 200
T = timebasisshiftedlagrange(Δt, Nt, 1)
U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U

duration = 20 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)

direction, polarisation = ẑ, x̂
E = BEAST.planewave(polarisation, direction, gaussian, 1.0)

SL = TDMaxwell3D.singlelayer(speedoflight=1.0)

@hilbertspace j
@hilbertspace j′
efie_nodot = @discretise SL[j′,j] == E[j′] j∈V j′∈W
# error()

xefie_nodot = solve(efie_nodot)

Xefie, Δω, ω0 = fouriertransform(xefie_nodot, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω.-1.0))

ω1 = ω[i1]
ue = Xefie[:,i1] / fouriertransform(gaussian)(ω1)
