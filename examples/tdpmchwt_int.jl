# using Pkg
# Pkg.activate(@__DIR__)
# Pkg.instantiate()

using CompScienceMeshes, BEAST
using LinearAlgebra

# Γ = meshcuboid(1.0, 1.0, 1.0, 0.25)
# h = 0.3
fn = joinpath(dirname(pathof(CompScienceMeshes)),"geos/torus.geo")
Γ = CompScienceMeshes.meshgeo(fn; dim=2, h=1.0)

Γ = meshsphere(radius=1.0, h=0.25)
X = raviartthomas(Γ)
Y = buffachristiansen(Γ, ibscaled=false)

sol = 1.0
T = TDMaxwell3D.singlelayer(speedoflight=sol, numdiffs=0)
K = TDMaxwell3D.doublelayer(speedoflight=sol, numdiffs=0)

Δt, Nt = 0.1, 300
δ = timebasisdelta(Δt, Nt)
T0 = timebasiscxd0(Δt, Nt)
T1 = timebasisshiftedlagrange(Δt,Nt,1)
T2 = timebasisshiftedlagrange(Δt, Nt, 2)
T3 = timebasisshiftedlagrange(Δt, Nt, 3)

duration = 40 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)

direction, polarisation = ẑ, x̂
E = BEAST.planewave(polarisation, direction, gaussian, 1.0)
H = direction × E

@hilbertspace j m
@hilbertspace k l

# pmchwt = @discretise(
#     2.0T[k,j] + 2.0K[k,m] -
#     2.0K[l,j] + 2.0T[l,m] == H[k] + E[l],
#     j∈X⊗T2, m∈Y⊗T2, k∈X⊗δ, l∈Y⊗δ)

pmchwt = @discretise(
    2.0T[k,j] + 2.0K[k,m] -
    2.0K[l,j] + 2.0T[l,m] == H[k] + E[l],
    j∈X⊗T1, m∈X⊗T1, k∈X⊗δ, l∈X⊗δ)

# Z = BEAST.td_assemble(pmchwt.equation.lhs, pmchwt.test_space_dict, pmchwt.trial_space_dict);
# error()


function cones(mesh, arrows; sizeref=2)
    centers = [cartesian(CompScienceMeshes.center(chart(mesh,cell))) for cell in mesh]
    x = getindex.(centers,1)
    y = getindex.(centers,2)
    z = getindex.(centers,3)
    u = getindex.(arrows,1)
    v = getindex.(arrows,2)
    w = getindex.(arrows,3)
    Plotly.cone(x=x,y=y,z=z,u=u,v=v,w=w,sizemode="absolute", sizeref=sizeref)
end



u = solve(pmchwt)
error()

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

# Study the pmchwt static nullspace
pmchwt = @discretise(
    2.0T[k,j] + 2.0K[k,m] -
    2.0K[l,j] + 2.0T[l,m] == H[k] + E[l],
    j∈X⊗T2, m∈X⊗T2, k∈X⊗δ, l∈X⊗δ)


# Z = BEAST.td_assemble(pmchwt.equation.lhs, pmchwt.test_space_dict, pmchwt.trial_space_dict)
# function cast(Z)
#     # kmax = maximum(length.(Z))
#     kmax = size(Z,3)
#     Q = zeros(eltype(eltype(Z)), size(Z)..., kmax)
#     for m in axes(Q,1)
#         for n in axes(Q,2)
#             for k in eachindex(Z[m,n])
#                 Q[m,n,k] = Z[m,n][k]
#             end
#         end
#     end
#     return Q
# end
# Q = cast(Z)

# C = companion(Q)
# w = eigvals(C)

# using Plots
# plotly()
# plot(exp.(im*range(0,2pi,length=200)))
# scatter!(w)