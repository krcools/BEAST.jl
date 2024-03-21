using DrWatson
quickactivate("nonlinearmodeling")
using CompScienceMeshes, BEAST

Γ = meshsphere(1.0, 0.2)
X = BEAST.raviartthomas(Γ)
c = 1.0 #speed of light
Δt, Nt = 0.2, 200
V = BEAST.FiniteDiffTimeStep(X, Δt, Nt)

duration = 2 * 20 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)

direction, polarisation = ẑ, x̂
E = planewave(polarisation, direction, derive(gaussian), 1.0)
LaplaceEFIO(s::T) where {T} = MWSingleLayer3D(s/c, s*s/c, T(c))
kmax = 30
rho = 1.0001
method = BEAST.BE(Δt)
T = BEAST.FiniteDiffConvolutionQuadrature(LaplaceEFIO, method, Δt, kmax, rho)

@hilbertspace j
@hilbertspace j′
tdefie = @discretise T[j′,j] == -1E[j′]   j∈V  j′∈V
xefie_cq = solve(tdefie)
