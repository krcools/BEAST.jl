using CompScienceMeshes, BEAST
o, x, y, z = euclidianbasis(3)

# D, Δx = 1.0, 0.35
# Γ = meshsphere(D, Δx)

Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
X = raviartthomas(Γ)
#Δt, Nt = 0.08, 400
Δt, Nt = 0.6, 200
T = timebasisc0d1(Δt, Nt)
U = timebasiscxd0(Δt, Nt)

V = X ⊗ T
W = X ⊗ U

width, delay, scaling = 8.0, 12.0, 1.0
gaussian = creategaussian(width, delay, scaling)

direction, polarisation = ẑ, x̂
E = BEAST.planewave(polarisation, direction, derive(gaussian), 1.0)



@hilbertspace j; @hilbertspace j′
T = MWSingleLayerTDIO(1.0,-1.0,-1.0,2,0)
tdefie = @discretise T[j′,j] == -1E[j′]   j∈V  j′∈W
xefie = solve(tdefie)

Xefie, Δω, ω0 = fouriertransform(xefie, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω.-1.0))

ω1 = ω[i1]
ue = Xefie[:,i1]/ fouriertransform(gaussian)(ω1)
