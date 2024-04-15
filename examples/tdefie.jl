using CompScienceMeshes, BEAST, LinearAlgebra
Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
Γ = meshsphere(radius=1.0, h=0.25)

X = raviartthomas(Γ)

Δt = 0.1
Nt = 200
T = timebasisshiftedlagrange(Δt, Nt, 3)
U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U

duration = 2 * 20 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)
direction, polarisation = ẑ, x̂
E = planewave(polarisation, direction, derive(gaussian), 1.0)

@hilbertspace j
@hilbertspace j′

SL = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)
# BEAST.@defaultquadstrat (SL, W, V) BEAST.OuterNumInnerAnalyticQStrat(7)

tdefie = @discretise SL[j′,j] == -1.0E[j′]   j∈V  j′∈W
xefie = BEAST.motsolve(tdefie)

import Plots
Plots.plot(xefie[1,:])

import Plotly
fcr, geo = facecurrents(xefie[:,125], X)
Plotly.plot(patch(geo, norm.(fcr)))






Xefie, Δω, ω0 = fouriertransform(xefie, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω.-1.0))

ω1 = ω[i1]
ue = Xefie[:,i1] / fouriertransform(gaussian)(ω1)

fcr, geo = facecurrents(ue, X)
Plotly.plot(patch(geo, norm.(fcr)))