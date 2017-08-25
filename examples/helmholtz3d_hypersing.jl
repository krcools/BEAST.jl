using CompScienceMeshes, BEAST

o, x, y, z = euclidianbasis(3)

Γ = meshsphere(1.0, 0.11)
X = lagrangec0d1(Γ)
@show numfunctions(X)

κ = 1.0; γ = im*κ
a = -Helmholtz3D.hypersingular(gamma=γ)
b = Helmholtz3D.doublelayer(gamma=γ) - 0.5Identity()

uⁱ = Helmholtz3D.planewave(wavenumber=κ, direction=z)
f = strace(uⁱ,Γ)
g = ∂n(uⁱ)


@hilbertspace u
@hilbertspace v
eq1 = @discretise a[v,u] == g[v] u∈X v∈X
eq2 = @discretise b[v,u] == f[v] u∈X v∈X

x1 = solve(eq1)
x2 = solve(eq2)

fcr1, geo1 = facecurrents(x1, X)
fcr2, geo2 = facecurrents(x2, X)

include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))
p1 = patch(geo1, real.(norm.(fcr1)))
p2 = patch(geo2, real.(norm.(fcr2)))
