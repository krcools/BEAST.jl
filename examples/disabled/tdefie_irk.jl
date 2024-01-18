using CompScienceMeshes, BEAST, StaticArrays, LinearAlgebra, Plots

Γ = meshsphere(radius=1.0, h=0.35)
X = raviartthomas(Γ)
sol = 1.0
Δt, Nt = 10.0, 200

(A, b, c) = butcher_tableau_radau_2stages()
V = StagedTimeStep(X, c, Δt, Nt)
LaplaceEFIO(s::T) where {T}= MWSingleLayer3D(-s/sol, s*s/sol, T(sol))

duration, delay, amplitude = 2 * 20 * Δt, 2 * 30 * Δt, 1.0
gaussian = creategaussian(duration, delay, amplitude)

direction, polarisation = ẑ , x̂
E = planewave(polarisation, direction, BEAST.derive(gaussian), sol)
T = RungeKuttaConvolutionQuadrature(LaplaceEFIO, A, b, Δt, 10, 1.001);

@hilbertspace j
@hilbertspace j′
tdefie_irk = @discretise T[j′,j] == -1E[j′]   j∈V  j′∈V
xefie_irk = solve(tdefie_irk)

j = xefie_irk[1:2:end,:]

import Plots
Plots.plot(j[1,:])

import Plotly
fcr, geo = facecurrents(j[:,1], X)
Plotly.plot(patch(geo, norm.(fcr)))

Xefie_irk, Δω, ω0 = fouriertransform(j, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω.-1.0))
ω1 = ω[i1]
ue = Xefie_irk[:, i1] / fouriertransform(gaussian)(ω1)
fcr, geo = facecurrents(ue, X)
Plotly.plot(patch(geo, norm.(fcr)))
