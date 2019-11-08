# This script solved the TD-EFIE and TD-MFIE 'at the same time'. It verifies
# the code dealing with systems of time domain equations in the special case
# of a block diagonal system.
using CompScienceMeshes, BEAST
Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))

X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

Δt, Nt = 0.6, 200
T2 = timebasisshiftedlagrange(Δt, Nt, 2)
δ = timebasisdelta(Δt, Nt)
T3 = timebasisshiftedlagrange(Δt, Nt, 3)

# V = X ⊗ T
# W = Y ⊗ δ

duration = 20 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)

direction, polarisation = ẑ, x̂
E = BEAST.planewave(polarisation, direction, gaussian, 1.0)
dE = BEAST.planewave(polarisation, direction, derive(gaussian), 1.0)
H = direction × E

dSL = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)
DL = TDMaxwell3D.doublelayer(speedoflight=1.0)
I = Identity()
N = NCross()

@hilbertspace k l
@hilbertspace j m
emfie = @discretise(
    (0.5(N⊗I) + 1.0DL)[k,j] + dSL[l,m] == -1.0H[k] - dE[l],
    k∈Y⊗δ, l∈X⊗δ, j∈X⊗T2, m∈X⊗T3)
xemfie = solve(emfie)

# Xmfie, Δω, ω0 = fouriertransform(xmfie, Δt, 0.0, 2)
# ω = collect(ω0 .+ (0:Nt-1)*Δω)
# _, i1 = findmin(abs.(ω .- 1.0))
#
# ω1 = ω[i1]
# um = Xmfie[:,i1] / fouriertransform(gaussian)(ω1)
