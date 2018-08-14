using CompScienceMeshes
using BEAST

G = meshsphere(1.0, 0.25)

c = 1.0
S = BEAST.HH3DSingleLayerTDBIO(c)

width, delay, scaling = 8.0, 12.0, 1.0
gaussian = creategaussian(width, delay, scaling)
e = BEAST.planewave(point(0,0,1), c, gaussian)

X = lagrangecxd0(G)

Δt, Nt = 0.08, 300
#T = timebasisshiftedlagrange(Δt, Nt, 0)
#T = timebasiscxd0(Δt, Nt)
T = timebasisc0d1(Δt, Nt)
U = timebasisdelta(Δt, Nt)

V = X ⊗ T # trial space
W = X ⊗ U # test space

b = assemble(e, W)
Z = assemble(S, W, V)

W = inv(Z[:,:,1])
u = marchonintime(W,Z,-b,Nt)

U, Δω, ω0 = fouriertransform(u, Δt, 0.0, 2)
ω = collect(ω0 + (0:Nt-1)*Δω)
_, i1 = findmin(abs(ω-1.0)); ω1 = ω[i1]

U1 = U[:,i1]
fgaussian = fouriertransform(gaussian)
U1 /= fgaussian(ω1)
Fcr, geo = facecurrents(U1, X)

#A = [real(norm(f)) for f in Fcr]
include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))
p = patch(geo, real.(norm.(fcr)))
#PlotlyJS.plot([p])
