using CompScienceMeshes
using BEAST
using LinearAlgebra
import Plots

κ, η = 1.0, 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

Φ, Θ = [0.0], range(0,stop=π,length=100)
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]
xs, y, zs = range(-4,stop=6,length=200), 1.0, range(-5,stop=5,length=200)
gridpoints = [point(x,y,z) for z in zs, x in xs]

Γ = CompScienceMeshes.meshrectangle(2.0, 2.0, 0.025; structured=:quadrilateral)
X = raviartthomas(Γ)

@hilbertspace j
@hilbertspace k
efie = @discretise t[k,j]==e[k]  j∈X k∈X
u, ch = BEAST.gmres_ch(efie; restart=1500)

ffd = potential(MWFarField3D(wavenumber=κ), pts, u, X)
nfd = potential(MWSingleLayerField3D(wavenumber = κ), gridpoints, u, X)
nfd .-= E.(gridpoints)

p1 = Plots.scatter(Θ, real.(norm.(ffd)))
p2 = Plots.heatmap(clamp.(abs.(getindex.(nfd,1)), 0, 2.0), clims=(0,2), colormap=:viridis)

Γ = CompScienceMeshes.meshrectangle(2.0, 2.0, 0.025)
X = raviartthomas(Γ)

@hilbertspace j
@hilbertspace k
efie = @discretise t[k,j]==e[k]  j∈X k∈X
u, ch = BEAST.gmres_ch(efie; restart=1500)

ffd = potential(MWFarField3D(wavenumber=κ), pts, u, X)
nfd = potential(MWSingleLayerField3D(wavenumber = κ), gridpoints, u, X)
nfd .-= E.(gridpoints)

p3 = Plots.scatter(Θ, real.(norm.(ffd)))
p4 = Plots.heatmap(clamp.(abs.(getindex.(nfd,1)), 0, 2.0), clims=(0,2), colormap=:viridis)

Plots.plot(p2,p4,layout=(1,2))

