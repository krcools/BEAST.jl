using CompScienceMeshes, BEAST
using Base.Test

o, x, y, z = euclidianbasis(3)

# Γ = meshsphere(1.0, 0.11)
Γ = readmesh("/Users/Benjamin/Documents/sphere.in")
X = lagrangecxd0(Γ)
@show numfunctions(X)

κ = 1.0; γ = im*κ
a = Helmholtz3D.singlelayer(gamma=γ)
b = Helmholtz3D.doublelayer_transposed(gamma=γ) +0.5Identity()

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

include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))
p1 = patch(geo1, real.(norm.(fcr1)))
p2 = patch(geo1, real.(norm.(fcr2)))

## test the results
Z = assemble(a,X,X);
m, n = 1, numfunctions(X)
chm, chn = chart(Γ,cells(Γ)[m]), chart(Γ,cells(Γ)[n])
ctm, ctn = center(chm), center(chn)
R = norm(cartesian(ctm)-cartesian(ctn))
G = exp(-im*κ*R)/(4π*R)
Wmn = volume(chm) * volume(chn) * G
@show abs(Wmn-Z[m,n]) / abs(Z[m,n])
@test abs(Wmn-Z[m,n]) / abs(Z[m,n]) < 2.0e-3

r = assemble(f,X)
m = 1
chm = chart(Γ,cells(Γ)[m])
ctm = center(chm)
sm = volume(chm) * f(ctm)
r[m]
@show abs(sm - r[m]) / abs(r[m])
@test abs(sm - r[m]) / abs(r[m]) < 1e-3
