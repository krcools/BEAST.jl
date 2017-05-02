using BEAST
using CompScienceMeshes
#using LinearForms
include(Pkg.dir("BEAST","src","lusolver.jl"))

o, x, y, z = euclidianbasis(3)
n = BEAST.n

Γ = meshsphere(1.0, 0.2)
RT = raviartthomas(Γ)

κ = 1.0
t = MWSingleLayer3D(im*κ)
E = planewavemw3d(direction=z, polarization=x, wavenumber=κ)
e = (n × E) × n

j, = hilbertspace(:j)
k, = hilbertspace(:k)

EFIE = @varform t[k,j] == e[k]
efie = @discretise EFIE j∈RT k∈RT

u = solve(efie)
println("Solution computed")

fcr, geo = facecurrents(u, RT)
println("Face currents computed")

include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))
A = real.(norm.(fcr))
p = patch(geo, A)
display(PlotlyJS.plot(p))

isdefined(:nearfar) || (nearfar = false;)
if nearfar
    pts = 2 * [point(cos(p)*sin(t), sin(p)*sin(t), cos(t)) for t in linspace(0,π,100) for p in [0.0]]
    ffd = potential(MWFarField3D(im*κ), pts, u, RT)
    println("far field computed")

    grid = [point(x,0,z) for x in linspace(-2,2,40), z in linspace(-4,4,40)]
    nfd = potential(MWSingleLayerField3D(κ), grid, u, RT)
    nfd = reshape(nfd, (40,40))
    nfd .-= E.(grid)
    println("near field computed.")

    A = [real(f[1]) for f in nfd]
    A = clamp(A,0.0, 2.5)

    using Plots
    heatmap(A)
end
