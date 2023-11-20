using CompScienceMeshes, BEAST, LinearAlgebra
Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
Γ = meshsphere(radius=1.0, h=0.5)

X = lagrangecxd0(Γ)
numfunctions(X)
Δt = 0.05
Nt = 20
T = timebasisshiftedlagrange(Δt, Nt, 0)
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

SL = TDAcustic3D.acusticsinglelayer(speedofsound=1.0, numdiffs=1)
# BEAST.@defaultquadstrat (SL, W, V) BEAST.OuterNumInnerAnalyticQStrat(7)

tdacusticsl = @discretise SL[j′,j] == -1.0E[j′]   j∈V  j′∈W
xacusticsl = solve(tdacusticsl)
#corregere da qui 
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