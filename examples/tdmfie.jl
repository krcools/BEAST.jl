using CompScienceMeshes, BEAST
Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
Γ = meshsphere(radius=1.0, h=0.25)

X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

Δt, Nt = 0.1, 200
T = timebasisshiftedlagrange(Δt, Nt, 2)
δ = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = Y ⊗ δ

duration = 2 * 20 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)

direction, polarisation = ẑ, x̂
E = BEAST.planewave(polarisation, direction, gaussian, 1.0)
H = direction × E

DL = TDMaxwell3D.doublelayer(speedoflight=1.0)
I = Identity()
N = NCross()

@hilbertspace k
@hilbertspace j
mfie = @discretise ((0.5*N)⊗I)[k,j] + 1.0DL[k,j] == -1.0H[k] k∈W j∈V

xmfie = BEAST.motsolve(mfie)

Xmfie, Δω, ω0 = fouriertransform(xmfie, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω .- 1.0))

ω1 = ω[i1]
um = Xmfie[:,i1] / fouriertransform(gaussian)(ω1)

import Plotly
fcr, geo = facecurrents(um, X)
Plotly.plot(patch(geo, norm.(fcr)))
