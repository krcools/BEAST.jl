using CompScienceMeshes, BEAST, StaticArrays, LinearAlgebra #, Plots

Γ = meshsphere(radius=1.0, h=0.45)
X = raviartthomas(Γ)
sol = 1.0
Δt, Nt = 10.0, 200

(A, b, c) = butcher_tableau_radau_3stages()
T = StagedTimeStep(Δt, Nt, c, A, b, 10, 1.0001)
V = X ⊗ T

duration = 2 * 20 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)

direction, polarisation = ẑ , x̂
E = planewave(polarisation, direction, derive(gaussian), sol)

T = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)

@hilbertspace j
@hilbertspace j′
tdefie_irk = @discretise T[j′,j] == -1E[j′]   j∈V  j′∈V
xefie_irk = solve(tdefie_irk)

import PlotlyJS
PlotlyJS.plot(xefie_irk[1,:])



