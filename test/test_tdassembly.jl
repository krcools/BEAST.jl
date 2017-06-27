using BEAST
using CompScienceMeshes
using WiltonInts84
using Base.Test


p1 = point(0,0,0)
p2 = point(1,0,0)
p3 = point(0,1,0)

q1 = point(2,2,0)
q2 = point(2,1,0)
q3 = point(1,2,0)

τ = [p1,p2,p3]
σ = [q1,q2,q3]
m, M = BEAST.minmaxdist(τ, σ)

@show m, sqrt(2.0)
@show M, sqrt(8.0)

@test m ≈ sqrt(2.0)
@test M ≈ sqrt(8.0)

ΔR = sqrt(2)/2
r = BEAST.rings(τ,σ,ΔR)
@show r

cm = point(0.5,0.5,0.0)
cM = point(0.0,0.0,0.0)
S,V,G = WiltonInts84.wiltonints(q1,q2,q3,cm,(r[1]-2)*ΔR,(r[1]-1)*ΔR,Val{0})
@show S[3]

S,V,G = WiltonInts84.wiltonints(q1,q2,q3,cM,(r[end]-1)*ΔR,(r[end])*ΔR,Val{0})
@show S[3]

h = 0.05
Γ = meshrectangle(10.0, h, h, 3)
X = raviartthomas(Γ)

timestep = 1.0
numsteps = 10
S = timebasisc0d1(timestep, numsteps)
S = BEAST.timebasisspline2(timestep, numsteps)

c = 1.0
Ts = MWSingleLayerTDIO(c,1.0,0.0,0,0)
Th = MWSingleLayerTDIO(c,0.0,1.0,0,0)

Zs = BEAST.assemble(Ts,X,X,S)
Zh = BEAST.assemble(Th,X,X,S)

## Test the allcatestorage routine and in particular its type stability
h = 0.05
Γ = meshrectangle(10.0, h, h, 3)
X = raviartthomas(Γ)

timestep = 1.0
numsteps = 10
S = timebasisc0d1(timestep, numsteps)
S = BEAST.timebasisspline2(timestep, numsteps)

c = 1.0
Ts = MWSingleLayerTDIO(c,1.0,0.0,0,0)
Th = MWSingleLayerTDIO(c,0.0,1.0,0,0)

U = timebasiscxd0(Δt, Nt)
T = timebasisc0d1(Δt, Nt)

W = X ⊗ U
V = X ⊗ T
@code_warntype BEAST.allocatestorage(Ts,W,V)
