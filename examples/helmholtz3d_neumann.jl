using CompScienceMeshes, BEAST
using LinearAlgebra, Pkg

Pkg.activate(@__DIR__)
# Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
X = lagrangec0d1(Γ)
@show numfunctions(X)

κ = 1.0; γ = im*κ
a = -1Helmholtz3D.hypersingular(gamma=γ)
b = Helmholtz3D.doublelayer(gamma=γ) - 0.5Identity()

uⁱ = Helmholtz3D.planewave(wavenumber=κ, direction=ẑ)
f = strace(uⁱ,Γ)
g = ∂n(uⁱ)

@hilbertspace u
@hilbertspace v
eq1 = @discretise a[v,u] == g[v] u∈X v∈X
eq2 = @discretise b[v,u] == f[v] u∈X v∈X

x1 = solve(eq1)
x2 = solve(eq2)

@assert norm(x1-x2)/norm(x1+x2) <  0.5e-2

fcr1, geo1 = facecurrents(x1, X)
fcr2, geo2 = facecurrents(x2, X)

using Plots
plot(title="Comparse 1st and 2nd kind eqs.")
plot!(norm.(fcr1),c=:blue,label="1st")
scatter!(norm.(fcr2),c=:red,label="2nd")

import Plotly
Plotly.plot(patch(Γ, norm.(fcr1)))
Plotly.plot(patch(Γ, norm.(fcr2)))