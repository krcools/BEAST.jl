using CompScienceMeshes
using BEAST

o, x, y, z = euclidianbasis(3)

D, Δx = 1.0, 0.25
Γ = meshsphere(D, Δx)
X = raviartthomas(Γ)


Δt, Nt = 0.08, 400
T = timebasisc0d1(Δt, Nt)
U = timebasiscxd0(Δt, Nt)

V = X ⊗ T
W = X ⊗ U

k = z
p = x
width = 8.0
delay = 12.0
scaling = 1.0
gaussian = creategaussian(width, delay, scaling)
fgaussian = fouriertransform(gaussian)
dgaussian = derive(gaussian)
E = planewave(p, k, derive(gaussian), 1.0)

T = MWSingleLayerTDIO(1.0,-1.0,-1.0,2,0)





B = assemble(E, W)
Z = assemble(T, W, V)

Z0 = Z[:,:,1]
W0 = inv(Z0)
tefie = 0.0 : Δt : (Nt-1)*Δt
xefie = marchonintime(W0,Z,-B,Nt)

Xefie, Δω, ω0 = fouriertransform(xefie, Δt, 0.0, 2)
ω = collect(ω0 + (0:Nt-1)*Δω)
_, i1 = findmin(abs(ω-1.0))
ω1 = ω[i1]

ue = Xefie[:,i1]
ue /= fgaussian(ω1)
fcre, geo = facecurrents(ue, X)

include(Pkg.dir("CompScienceMeshes","examples","matlab_patches.jl"))
p = patch(geo, real.(norm.(fcre)))
PlotlyJS.plot([p])
