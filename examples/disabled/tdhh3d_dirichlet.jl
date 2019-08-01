using CompScienceMeshes
using BEAST

G = meshsphere(1.0, 0.25)

c = 1.0
# S = BEAST.HH3DSingleLayerTDBIO(c)
S = BEAST.HH3DHyperSingularTDBIO(speed_of_light=c, numdiffs=1)

width, delay, scaling = 8.0, 12.0, 1.0
gaussian = creategaussian(width, delay, scaling)
e = BEAST.planewave(point(0,0,1), c, gaussian)

# X = lagrangecxd0(G)
X = lagrangec0d1(G)

Δt, Nt = 0.08, 300
T = timebasisc0d1(Δt, Nt)
U = timebasiscxd0(Δt, Nt)

V = X ⊗ T # trial space
W = X ⊗ U # test space

nbd = center(chart(G, first(cells(G))))
refs = refspace(X)
vals = refs(nbd)

Z = assemble(S, W, V)
h = dot(n,BEAST.gradient(e))
b = assemble(h, W)

iZ1 = inv(Z[:,:,1])
u = marchonintime(iZ1,Z,-b,Nt)

U, Δω, ω0 = fouriertransform(u, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω.-1.0)); ω1 = ω[i1]

U1 = U[:,i1]
fgaussian = fouriertransform(gaussian)
U1 /= fgaussian(ω1)
Fcr, geo = facecurrents(U1, X)

include(joinpath(dirname(pathof(CompScienceMeshes)),"..","examples","plotlyjs_patches.jl"))

using LinearAlgebra
p = patch(geo, real.(norm.(Fcr)))
#PlotlyJS.plot([p])
