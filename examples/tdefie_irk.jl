include("../src/BEAST.jl")
using CompScienceMeshes, BEAST, StaticArrays

o, x, y, z = euclidianbasis(3)

sol = 5.0;
#Δt, Nt = 0.08/sol,400
Δt, Nt = 100.0/sol,200

D, Δx = 1.0, 0.45
Γ = meshsphere(D, Δx)
X = raviartthomas(Γ)

(A, b, c) = butcher_tableau_radau_2stages();
V = StagedTimeStep(X, c, Δt, Nt);

duration, delay, amplitude = 2000.0/sol, 2000.0/sol, 1.0
gaussian = creategaussian(duration, delay, duration)

direction, polarisation = z, x
E = planewave(polarisation, direction, derive(gaussian), sol)

LaplaceEFIO{T}(s::T) = MWSingleLayer3D(-s/sol, s*s/sol, T(sol));
T = RungeKuttaConvolutionQuadrature(LaplaceEFIO, A, b, Δt, 15, 1.0001);

@hilbertspace j
@hilbertspace j′
tdefie = @discretise T[j′,j] == -1E[j′]   j∈V  j′∈V
xefie_irk = solve(tdefie)

using PyPlot
semilogy(abs.(xefie_irk[2,:]))
