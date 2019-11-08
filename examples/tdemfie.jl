# This script solved the TD-EFIE and TD-MFIE 'at the same time'. It verifies
# the code dealing with systems of time domain equations in the special case
# of a block diagonal system.
using CompScienceMeshes, BEAST
Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))

X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

Δt, Nt = 0.6, 200
δ = timebasisdelta(Δt, Nt)
T0 = timebasiscxd0(Δt, Nt)
T1 = timebasisshiftedlagrange(Δt,Nt,1)
T2 = timebasisshiftedlagrange(Δt, Nt, 2)
T3 = timebasisshiftedlagrange(Δt, Nt, 3)

duration = 20 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)

direction, polarisation = ẑ, x̂
E = BEAST.planewave(polarisation, direction, gaussian, 1.0)
H = direction × E

SL = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=0)
DL = TDMaxwell3D.doublelayer(speedoflight=1.0)
I = Identity()
N = NCross()

@hilbertspace k l
@hilbertspace j m
emfie = @discretise(
    (0.5(N⊗I) + 1.0DL)[k,j] + SL[l,m] == -1.0H[k] - E[l],
    k∈Y⊗δ, l∈X⊗δ, j∈X⊗T1, m∈X⊗T1)
xemfie = solve(emfie)
