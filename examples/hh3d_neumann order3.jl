using CompScienceMeshes, BEAST
using LinearAlgebra, Pkg

# Pkg.activate(@__DIR__)
Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
# Γ = meshrectangle(1.0, 1.0, 0.05, 3)
X = BEAST.lagrangec0(Γ; order=2)
@show numfunctions(X)

κ = 2π; γ = im*κ
a = Helmholtz3D.hypersingular(gamma=γ)
# b = Helmholtz3D.doublelayer(gamma=γ) - 0.5Identity()

uⁱ = Helmholtz3D.planewave(wavenumber=κ, direction=ẑ)
# f = strace(uⁱ,Γ)
g = ∂n(uⁱ)

BEAST.@defaultquadstrat (a,X,X) BEAST.DoubleNumSauterQstrat(7,8,6,6,6,6)
BEAST.@defaultquadstrat (g,X) BEAST.SingleNumQStrat(12)

@hilbertspace u
@hilbertspace v

A = assemble(a[v,u], X, X)
b = assemble(g[v], X)
x1 = AbstractMatrix(A) \ b

# eq1 = @discretise a[v,u] == g[v] u∈X v∈X
# eq2 = @discretise b[v,u] == f[v] u∈X v∈X

# x1 = gmres(eq1)
# x2 = gmres(eq2)

# @assert norm(x1-x2)/norm(x1+x2) <  0.5e-2

fcr1, geo1 = facecurrents(x1, X)
# fcr2, geo2 = facecurrents(x2, X)

using Plots
Plots.plot(title="Comparse 1st and 2nd kind eqs.")
Plots.plot!(norm.(fcr1),c=:blue,label="1st")
# Plots.scatter!(norm.(fcr2),c=:red,label="2nd")

import Plotly
plt1 = Plotly.plot(patch(Γ, norm.(fcr1)))
# plt2 = Plotly.plot(patch(Γ, norm.(fcr2)))
# display([plt1 plt2])

# ys = range(-2,2,length=200)
# zs = range(-2,2,length=200)
# pts = [point(0.5, y, z) for y in ys, z in zs]

# nf = BEAST.HH3DDoubleLayerNear(wavenumber=κ)
# near = BEAST.potential(nf, pts, x1, X, type=ComplexF64)
# inc = uⁱ.(pts)
# tot = near + inc
