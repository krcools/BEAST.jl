using CompScienceMeshes, BEAST

o, x, y, z = euclidianbasis(3)

# Γ = meshsphere(1.0, 0.11)
Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
# Γ = readmesh("/Users/Benjamin/Documents/sphere.in")
# Γ = readmesh(joinpath(@__DIR__,"sphere_subd1.in"))
X = lagrangecxd0(Γ)
# X = subdsurface(Γ)
# X = raviartthomas(Γ)
@show numfunctions(X)

κ = 1.0; γ = im*κ
a = Helmholtz3D.singlelayer(wavenumber=κ)
b = Helmholtz3D.doublelayer_transposed(gamma=κ*im) +0.5Identity()

uⁱ = Helmholtz3D.planewave(wavenumber=κ, direction=z)
f = strace(uⁱ,Γ)
g = ∂n(uⁱ)

@hilbertspace u
@hilbertspace v

eq1 = @discretise a[v,u] == f[v] u∈X v∈X
eq2 = @discretise b[v,u] == g[v] u∈X v∈X

x1 = solve(eq1)
x2 = solve(eq2)

fcr1, geo1 = facecurrents(x1, X)
fcr2, geo2 = facecurrents(x2, X)

# include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))
# p1 = patch(geo1, real.(norm.(fcr1)))
# p2 = patch(geo2, real.(norm.(fcr2)))

using LinearAlgebra
using Test

## test the results
Z = assemble(a,X,X);
m1, m2 = 1, numfunctions(X)
chm, chn = chart(Γ,cells(Γ)[m1]), chart(Γ,cells(Γ)[m2])
ctm, ctn = center(chm), center(chn)
R = norm(cartesian(ctm)-cartesian(ctn))
G = exp(-im*κ*R)/(4π*R)
Wmn = volume(chm) * volume(chn) * G
@show abs(Wmn-Z[m1,m2]) / abs(Z[m1,m2])
@test abs(Wmn-Z[m1,m2]) / abs(Z[m1,m2]) < 2.0e-3

r = assemble(f,X)
m1 = 1
chm = chart(Γ,cells(Γ)[m1])
ctm = center(chm)
sm = volume(chm) * f(ctm)
r[m1]
@show abs(sm - r[m1]) / abs(r[m1])
@test abs(sm - r[m1]) / abs(r[m1]) < 1e-3
