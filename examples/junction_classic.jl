using CompScienceMeshes, BEAST
o, x, y, z = euclidianbasis(3)

width, height, h = 1.0, 0.5, 0.05
G1 = meshrectangle(width, height, h)
G2 = CompScienceMeshes.rotate(G1, 0.5π * [1,0,0])
G3 = CompScienceMeshes.rotate(G1, 1.0π * [1,0,0])
G = weld(G1, G2, G3)

X = raviartthomas(G)

κ = 1.0
E = planewavemw3d(direction=z, polarization=x, wavenumber=κ)
e = (n × E) × n
T = MWSingleLayer3D(-im*κ)

@hilbertspace j; @hilbertspace k
efie = @discretise T[k,j] == e[k]  j∈X k∈X
u = solve(efie)

Θ, Φ = linspace(0.0,π,7), 0.0
ffpoints = vec([point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ, ϕ in Φ])
ff = potential(MWFarField3D(κ*im), ffpoints, u, X)

include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))
fcr, geo = facecurrents(u, X)
p = patch(geo, real.(norm.(fcr)))
PlotlyJS.plot(p)
