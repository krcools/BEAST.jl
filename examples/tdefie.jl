using CompScienceMeshes, BEAST, LinearAlgebra
Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))

X = raviartthomas(Γ)

Δt, Nt = 0.3, 200
T = timebasisshiftedlagrange(Δt, Nt, 3)
U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U

duration = 20 * Δt * 2
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)
direction, polarisation = ẑ, x̂
E = planewave(polarisation, direction, derive(gaussian), 1.0)

@hilbertspace j
@hilbertspace j′

SL = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)
tdefie = @discretise SL[j′,j] == -1.0E[j′]   j∈V  j′∈W
xefie = solve(tdefie)

import Plots
Plots.plot(xefie[1,:])

import Plotly
fcr, geo = facecurrents(xefie[:,60], X)
Plotly.plot(patch(geo, norm.(fcr)))






Xefie, Δω, ω0 = fouriertransform(xefie, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω.-1.0))

ω1 = ω[i1]
ue = Xefie[:,i1] / fouriertransform(gaussian)(ω1)
