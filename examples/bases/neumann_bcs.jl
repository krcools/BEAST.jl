using BEAST
using CompScienceMeshes


Γ = meshrectangle(1.0, 1.0, 0.1, 3)
all_edges = skeleton(Γ, 1)
bnd_edges = boundary(Γ)
int_edges = setminus(all_edges, bnd_edges)

X = raviartthomas(Γ, all_edges)
Y = Ydi = BEAST.buffachristiansen(Γ, bnd_edges)

T = Maxwell3D.singlelayer(wavenumber=1.0)
Id = BEAST.NCross()

Txx = assemble(T, X, X)
Tyy = assemble(T, Y, Y)

Ixy = assemble(Id, X, Y)
IYX = BEAST.GMRES(Ixy)
IXY = BEAST.GMRES(Ixy')

using LinearAlgebra
@show cond(Txx)

𝗧 = Matrix(Txx)
𝗜 = Matrix(Ixy)
𝗜⁻¹ = inv(𝗜)
𝗣 = transpose(𝗜⁻¹) * Matrix(Tyy) * 𝗜⁻¹
𝗣𝗧 = 𝗣*𝗧

@show cond(𝗧)
@show cond(𝗜)
@show cond(𝗣)
@show cond(𝗣*𝗧)

σI = svdvals(𝗜)

using PlotlyJS
plot(
    scatter(y=σI),
    Layout(
        yaxis_type=:log))

