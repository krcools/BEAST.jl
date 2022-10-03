using CompScienceMeshes, BEAST
o, x, y, z = euclidianbasis(3)

sol = 5.0
#Δt, Nt = 0.08/sol,400
Δt, Nt = 0.120000123/sol,400

#D, Δx = 1.0, 0.25
D, Δx = 1.0, 0.1
#Γ = meshsphere(D, Δx)
Γ = meshrectangle(D,D,Δx,3)
#γ = boundary(Γ)
γ1 = meshsegment(D,D,3)
γ2 = CompScienceMeshes.translate(γ1, point(0,D,0))

#X = raviartthomas(Γ,weld(γ1,γ2))
X = rt_ports(Γ,(γ1,γ2))
X.fns = X.fns[1:281]
numfunctions(X)

T = timebasisshiftedlagrange(Δt, Nt, 3)
U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U

duration, delay, amplitude = 8.0/sol, 12.0/sol, 1.0
gaussian = creategaussian(duration, delay, duration)

direction, polarisation = z, x
E = planewave(polarisation, direction, derive(gaussian), sol)
T = MWSingleLayerTDIO(sol,-1/sol,-sol,2,0)

#Z = assemble(T,W,V);

@hilbertspace j
@hilbertspace j′
tdefie = @discretise T[j′,j] == -1E[j′]   j∈V  j′∈W
xefie = solve(tdefie)

using Plotly
# include(Pkg.dir("CompScienceMeshes","examples","matlab_patches.jl"))


Xefie, Δω, ω0 = fouriertransform(xefie, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω.-sol))

ω1 = ω[i1]
ue = Xefie[:,i1] / fouriertransform(gaussian)(ω1)

using LinearAlgebra
fcre, geo = facecurrents(ue, X)
t2 = patch(geo, real.(norm.(fcre)))
Plotly.plot(t2)


fcr_td, geo = facecurrents(xefie[:,end], X)
# fcr_ch, geo = facecurrents(xefie[:,end], divergence(X))

t3 = patch(Γ,real.(norm.(fcr_td)))
Plotly.plot(t3)
# jmatlab_quiver(Γ, fcr_td)
