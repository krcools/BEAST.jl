using CompScienceMeshes, BEAST
using LinearAlgebra

Γ = meshcuboid(1.0, 1.0, 1.0, 0.25)
# Γ = meshsphere(1.0, 0.3)
X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

sol = 1.0
T = TDMaxwell3D.singlelayer(speedoflight=sol, numdiffs=1)
K = TDMaxwell3D.doublelayer(speedoflight=sol, numdiffs=1)

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
E = BEAST.planewave(polarisation, direction, derive(gaussian), 1.0)
H = direction × E

@hilbertspace j m
@hilbertspace k l

pmchwt = @discretise(
    2.0T[k,j] + 2.0K[k,m] -
    2.0K[l,j] + 2.0T[l,m] == H[k] + E[l],
    j∈X⊗T2, m∈Y⊗T2, k∈X⊗δ, l∈Y⊗δ)

u = solve(pmchwt)

# nX = numfunctions(X)
# uj = u[1:nX]
# um = u[nX+1:end]

# Θ, Φ = range(0.0,stop=2π,length=100), 0.0
# ffpoints = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]
#
# # Don't forgt the far field comprises two contributions
# κ, η = 1.0, 1.0
# ffm = potential(MWFarField3D(κ*im), ffpoints, um, X)
# ffj = potential(MWFarField3D(κ*im), ffpoints, uj, X)
# ff = η*im*κ*ffj + im*κ*cross.(ffpoints, ffm)
#
# using Plots
# plot(xlabel="theta")
# plot!(Θ,norm.(ffm),label="far field")
#
# import PlotlyJS
# using LinearAlgebra
# fcrj, _ = facecurrents(uj,X)
# PlotlyJS.plot(patch(Γ, norm.(fcrj)))
