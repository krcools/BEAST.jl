using CompScienceMeshes, BEAST
o, x, y, z = euclidianbasis(3)

Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
X, Y = raviartthomas(Γ), buffachristiansen(Γ)

Δt ,Nt = 0.3, 200
T = timebasisshiftedlagrange(Δt, Nt, 2)
δ = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = Y ⊗ δ
duration = 20 * Δt * 2
delay = 1.5 * duration
amplitude = 1.0
gaussian = derive(creategaussian(duration, delay, amplitude))

direction, polarisation = ẑ, x̂
Ė = BEAST.planewave(polarisation, direction, gaussian, 1.0)
Ḣ = direction × Ė

@hilbertspace j
@hilbertspace k
K̇ = TDMaxwell3D.doublelayer(speedoflight=1.0, numdiffs=1)
Nİ = BEAST.TemporalDifferentiation(NCross()⊗Identity())

@hilbertspace k
@hilbertspace j
mfie_dot = @discretise (0.5*Nİ)[k,j] + K̇[k,j] == -1.0Ḣ[k] k∈W j∈V
xmfie_dot = BEAST.motsolve(mfie_dot)

# Xmfie, Δω, ω0 = fouriertransform(xmfie, Δt, 0.0, 2)
# ω = collect(ω0 .+ (0:Nt-1)*Δω)
# _, i1 = findmin(abs.(ω .- 1.0))
#
# ω1 = ω[i1]
# um = Xmfie[:,i1] / fouriertransform(gaussian)(ω1)
