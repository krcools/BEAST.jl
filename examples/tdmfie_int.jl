using CompScienceMeshes, BEAST
# Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
# Γ = meshsphere(1.0, 0.4)
@show length(Γ)

X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

Δt, Nt = 0.6, 200
T0 = timebasisshiftedlagrange(Δt, Nt, 0)
T1 = timebasisshiftedlagrange(Δt, Nt, 1)
iT0 = integrate(T0)
iT1 = integrate(T1)
δ = timebasisdelta(Δt, Nt)

# V = X ⊗ T
# iV = X ⊗ iT
# W = Y ⊗ δ

duration = 20 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)
int_gaussian = integrate(gaussian)

direction, polarisation = ẑ, x̂
E = BEAST.planewave(polarisation, direction, gaussian, 1.0)
iE = BEAST.planewave(polarisation, direction, int_gaussian, 1.0)
H = direction × E
iH = direction × iE

DL = TDMaxwell3D.doublelayer(speedoflight=1.0)
I = Identity()
N = NCross()
M = 0.5*N⊗I + 1.0*DL
iM = integrate(M)
NI = N ⊗ I
iNI = integrate(NI)

# error()

@hilbertspace k
@hilbertspace j
# mfie = @discretise (0.5(N⊗I) + 1.0DL)[k,j] == -1.0H[k] k∈W j∈V
mfie = @discretise iM[k,j] == -1.0iH[k] k∈(Y⊗δ) j∈(X⊗T0)

xmfie_int = solve(mfie)

Xmfie, Δω, ω0 = fouriertransform(xmfie_int, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω .- 1.0))

ω1 = ω[i1]
um = Xmfie[:,i1] / fouriertransform(gaussian)(ω1)
