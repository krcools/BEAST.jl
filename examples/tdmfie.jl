using CompScienceMeshes
using BEAST

x = point(1.0,0.0,0.0)
y = point(0.0,1.0,0.0)
z = point(0.0,0.0,1.0)

D, Δx = 1.0, 0.25
Γ = meshsphere(D, Δx)
X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

Δt, Nt = 0.1, 600
T = timebasisshiftedlagrange(Δt, Nt, 2)
δ = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = Y ⊗ δ

k = z
p = x
width = 8.0
delay = 12.0
scaling = 1.0
gaussian = creategaussian(width, delay, scaling)
fgaussian = fouriertransform(gaussian)
E = planewave(p, k, gaussian, 1.0)
H = k × E

K  = MWDoubleLayerTDIO(1.0, 1.0, 0)
I  = Identity()
N  = NCross()
NI = N ⊗ I
M = 0.5*(N⊗I) + 1.0*K

b = assemble(H, W)
Z = assemble(M, W, V)

Z0 = Z[:,:,1]
W0 = inv(Z0)
tmfie = 0.0 : Δt : (Nt-1)*Δt
xmfie = marchonintime(W0, Z, -b, Nt)

Xmfie, Δω, ω0 = fouriertransform(xmfie, Δt, 0.0, 2)
ω = collect(ω0 + (0:Nt-1)*Δω)
_, i1 = findmin(abs(ω-1.0))
ω1 = ω[i1]

um = Xmfie[:,i1]
um /= fgaussian(ω1)
fcrm, geo = facecurrents(um, X)

include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))
#A = [real(norm(f)) for f in fcrm]
p = patch(geo, real.(norm.(fcrm)))
PlotlyJS.plot([p])
