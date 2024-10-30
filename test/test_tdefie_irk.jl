using CompScienceMeshes
using BEAST
using StaticArrays
using LinearAlgebra
using Test

# Γ = meshsphere(radius=1.0, h=0.45)
Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../test/assets/sphere45.in"))
X = raviartthomas(Γ)
sol = 1.0 # Speed of light (for sake of simplicity, set to one)

# Note that certain choices of Δt, Nt, as well as the excitation
# can lead to late-time instabilities
Δt, Nt = 10.0, 200

(A, b, c) = butcher_tableau_radau_2stages()
T = StagedTimeStep(Δt, Nt, c, A, b, 5, 1.0001)
V = X ⊗ T

duration = 2 * 20 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)

direction, polarisation = ẑ , x̂
E = planewave(polarisation, direction, derive(gaussian), sol)

T = TDMaxwell3D.singlelayer(speedoflight=sol, numdiffs=1)

@hilbertspace j
@hilbertspace j′
tdefie_irk = @discretise T[j′,j] == -1E[j′]   j∈V  j′∈V
xefie_irk = solve(tdefie_irk)

# Set up Space-Time-Galerkin MOT for comparison

T = timebasisshiftedlagrange(Δt, Nt, 3)
U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U

SL = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)

tdefie = @discretise SL[j′,j] == -1.0E[j′]   j∈V  j′∈W
xefie = BEAST.motsolve(tdefie)

xefie_irk = xefie_irk[1:size(c,1):end,:]

diff_MOT_max = maximum((norm.(xefie_irk[:,1:end] - xefie[:,1:end])) ./ maximum(xefie[:,1:end]))

# @test diff_MOT_max ≈ 0.137252874891817 atol=1e-8
@test diff_MOT_max < 0.16 