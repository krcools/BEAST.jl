using CompScienceMeshes, BEAST
o, x, y, z = euclidianbasis(3)

sol = 5.0
#Δt, Nt = 0.08/sol,400
Δt, Nt = 0.12/sol,400

#D, Δx = 1.0, 0.25
D, Δx = 1.0, 0.35
Γ = meshsphere(D, Δx)
X = raviartthomas(Γ)

T = timebasisshiftedlagrange(Δt, Nt, 3)
U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U

duration, delay, amplitude = 8.0/sol, 12.0/sol, 1.0
gaussian = creategaussian(duration, delay, duration)

direction, polarisation = z, x
E = planewave(polarisation, direction, derive(gaussian), sol)
T = MWSingleLayerTDIO(sol,-1/sol,-sol,2,0)

Z = assemble(T,W,V);

# @hilbertspace j
# @hilbertspace j′
# tdefie = @discretise T[j′,j] == -1E[j′]   j∈V  j′∈W
# xefie = solve(tdefie)

# using PlotlyJS
# include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))
#
#
# Xefie, Δω, ω0 = fouriertransform(xefie, Δt, 0.0, 2)
# ω = collect(ω0 + (0:Nt-1)*Δω)
# _, i1 = findmin(abs.(ω-1.0*sol))
#
# ω1 = ω[i1]
# ue = Xefie[:,i1] / fouriertransform(gaussian)(ω1)
#
# fcre, geo = facecurrents(ue, X)
# t2 = patch(geo, real.(norm.(fcre)))
# PlotlyJS.plot(t2)
