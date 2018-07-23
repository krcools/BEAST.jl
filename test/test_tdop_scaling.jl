using CompScienceMeshes, BEAST, Test

fn = joinpath(dirname(@__FILE__),"assets","sphere35.in")

sol = 1.0
D, Δx = 1.0, 0.35
#Γ = meshsphere(D, Δx)
Γ = readmesh(fn)
X = raviartthomas(Γ)
Xm = subset(X,[1])
Xn = subset(X,[100])
Δt, Nt = 0.13/sol, 40
T = timebasisshiftedlagrange(Δt, Nt, 3)
U = timebasisdelta(Δt, Nt)

W = Xm ⊗ U
V = Xn ⊗ T

t = MWSingleLayerTDIO(sol,-1/sol,-sol,2,0)
Z1 = assemble(t, W, V)

## Run again, this time scales with another speedOfLight

sol = 3.0
D, Δx = 1.0, 0.35
#Γ = meshsphere(D, Δx)
Γ = readmesh(fn)
X = raviartthomas(Γ)
Xm = subset(X,[1])
Xn = subset(X,[100])
Δt, Nt = 0.13/sol, 40
T = timebasisshiftedlagrange(Δt, Nt, 3)
U = timebasisdelta(Δt, Nt)

W = Xm ⊗ U
V = Xn ⊗ T

t = MWSingleLayerTDIO(sol,-1/sol,-sol,2,0)
Z3 = assemble(t, W, V)

m = n = 1
@test findall(Z1[m,n,:] .!= 0) == findall(Z3[m,n,:] .!= 0)

I = findall(Z1[m,n,:] .!= 0)
@test all(sol*Z1[m,n,I] .≈ Z3[m,n,I])
