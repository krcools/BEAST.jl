using CompScienceMeshes, BEAST
using LinearAlgebra

Γ = meshcuboid(1.0, 1.0, 1.0, 0.15)
X = raviartthomas(Γ)

sol = 1.0
T = TDMaxwell3D.singlelayer(speedoflight=sol)
K = TDMaxwell3D.doublelayer(speedoflight=sol)

Δt, Nt = 0.6, 200
tbf = timebasisshiftedlagrange(Δt, Nt, 2)
ttf = timebasisdelta(Δt, Nt)

V = X ⊗ tbf
W = X ⊗ ttf

duration = 20 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)

direction, polarisation = ẑ, x̂
E = BEAST.planewave(polarisation, direction, gaussian, 1.0)
H = direction × E

@hilbertspace j m
@hilbertspace k l

pmchwt = @discretise(
    2.0T[k,j] + 2.0K[k,m] -
    2.0K[l,j] + 2.0T[l,m] == H[k] + E[l],
    j∈V, m∈V, k∈W, l∈W)

error("stop")

u = solve(pmchwt)

nX = numfunctions(X)
uj = u[1:nX]
um = u[nX+1:end]

Θ, Φ = range(0.0,stop=2π,length=100), 0.0
ffpoints = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]

# Don't forgt the far field comprises two contributions
κ, η = 1.0, 1.0
ffm = potential(MWFarField3D(κ*im), ffpoints, um, X)
ffj = potential(MWFarField3D(κ*im), ffpoints, uj, X)
ff = η*im*κ*ffj + im*κ*cross.(ffpoints, ffm)

using Plots
plot(xlabel="theta")
plot!(Θ,norm.(ffm),label="far field")

import PlotlyJS
using LinearAlgebra
fcrj, _ = facecurrents(uj,X)
PlotlyJS.plot(patch(Γ, norm.(fcrj)))
