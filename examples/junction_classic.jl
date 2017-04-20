using CompScienceMeshes
using BEAST

width, height, h = 1.0, 0.5, 0.1

G1 = meshrectangle(width, height, h)
G2 = rotate(G1, 0.5π * [1,0,0])
G3 = rotate(G1, 1.0π * [1,0,0])

G = weld(G1, G2, G3)
Test.@test numcells(G) == numcells(G1) + numcells(G2) + numcells(G3)
Test.@test numvertices(G) == 3 * numvertices(G1) - 2 * 11

E1 = skeleton(G1,1)
E = skeleton(G, 1)
P1 = cellpairs(G1,E1)
P = cellpairs(G,E)

X1 = raviartthomas(G1)
X = raviartthomas(G)
Test.@test numfunctions(X) == 3 * numfunctions(X1) + 20

κ = 1.0
x = point(1.0, 0.0, 0.0)
y = point(0.0, 1.0, 0.0)
z = point(0.0, 0.0, 1.0)
n = BEAST.n

E = PlaneWaveMW(z, x, κ)
e = (n × E) × n
b = assemble(e, X)

T = MWSingleLayer3D(-im*κ)
Txx = assemble(T, X, X)

u = Txx \ b

fcr, geo = facecurrents(u, X)
A = real.(nor.(fcr))

Θ, Φ = linspace(0.0,π,7), 0.0
ffpoints = vec([point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ, ϕ in Φ])
ff = potential(MWFarField3D(κ*im), ffpoints, u, X)

include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))
p = patch(geo, A)
PlotlyJS.plot(p)
extrema(A)
