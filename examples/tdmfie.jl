using CompScienceMeshes, BEAST
o, x, y, z = euclidianbasis(3)

Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
@show numcells(Γ)

X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

Δt, Nt = 0.6, 200
T = timebasisshiftedlagrange(Δt, Nt, 2)
δ = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = Y ⊗ δ

duration = 20 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)

direction, polarisation = ẑ, x̂
E = BEAST.planewave(polarisation, direction, gaussian, 1.0)
H = direction × E

@hilbertspace j; @hilbertspace m′
K = TDMaxwell3D.doublelayer(speedoflight=1.0)
I = Identity()
N = NCross()

M = 0.5*(N⊗I) + 1.0*K
Z_mfie = assemble(M, W, V, Val{:bandedstorage})
b_mfie = assemble(H, W)
xmfie = marchonintime(inv(Z_mfie[:,:,1]), Z_mfie, b_mfie, Nt)

Xmfie, Δω, ω0 = fouriertransform(xmfie, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω .- 1.0))

ω1 = ω[i1]
um = Xmfie[:,i1] / fouriertransform(gaussian)(ω1)
