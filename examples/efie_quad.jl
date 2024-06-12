using CompScienceMeshes
using BEAST

Γ = CompScienceMeshes.meshrectangle(2.0, 2.0, 0.1; structured=:quadrilateral)
edges = skeleton(Γ, 1)
edges_bnd = boundary(Γ)
pred = !in(edges_bnd)
edges_int = submesh(pred,  edges)
c = CompScienceMeshes.connectivity(edges_int, Γ, identity)
o = ones(length(edges_int))
X = raviartthomas(Γ, edges_int, c, o)
@show numfunctions(X)

κ, η = 1.0, 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

@hilbertspace j
@hilbertspace k
efie = @discretise t[k,j]==e[k]  j∈X k∈X
u, ch = BEAST.gmres_ch(efie; restart=1500)

Φ, Θ = [0.0], range(0,stop=π,length=100)
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]
ffd = potential(MWFarField3D(wavenumber=κ), pts, u, X)

xs, y, zs = range(-4,stop=6,length=200), 1.0, range(-5,stop=5,length=200)
gridpoints = [point(x,y,z) for z in zs, x in xs]
nfd = potential(MWSingleLayerField3D(wavenumber = κ), gridpoints, u, X)
# nfd = reshape(nfd, (nx,nz))
nfd .-= E.(gridpoints)

import Plots
using LinearAlgebra
p1 = Plots.scatter(Θ, real.(norm.(ffd)))
p2 = Plots.heatmap(-clamp.(real.(getindex.(nfd,1)), -3.0, 3.0), clims=(-3,3), colormap=:viridis)
p3 = Plots.contour(-clamp.(real.(getindex.(nfd,1)), -3.0, 3.0), clims=(-3,3))
Plots.plot(p1,p2,p3,layout=(3,1))

