using CompScienceMeshes, BEAST

width, height, h = 1.0, 0.5, 0.05
G1 = meshrectangle(width, height, h)
G2 = CompScienceMeshes.rotate(G1, 0.5π * x̂)
G3 = CompScienceMeshes.rotate(G1, 1.0π * x̂)
G = weld(G1, G2, G3)

X = raviartthomas(G)

κ = 1.0
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n
T = Maxwell3D.singlelayer(wavenumber=κ)

@hilbertspace j; @hilbertspace k
efie = @discretise T[k,j]==e[k]  j∈X k∈X
u = gmres(efie)

Θ, Φ = range(0.0,stop=π,length=100), 0.0
ffpoints = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]
farfield = potential(MWFarField3D(κ*im), ffpoints, u, X)

using Plots
using LinearAlgebra
plot(Θ,norm.(farfield))
