using CompScienceMeshes
using BEAST

Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
X = raviartthomas(Γ)
Y = buffachristiansen(Γ, ibscaled=true)

N = BEAST.NCross()

κ = 1.0
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

b = assemble(e, X)
G = assemble(N, X, Y)

x = G \ b
fcr, geo = facecurrents(x, Y)

using Plotly
using LinearAlgebra
Plotly.plot(patch(geo.mesh, norm.(fcr), (0.0, 1.0)))
