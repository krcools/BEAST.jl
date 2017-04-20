using CompScienceMeshes
using BEAST

using Plots


# TM scattering
h = 2π / 51
Γ = meshcircle(1.0, h)

X = lagrangecxd0(Γ)

κ = ω = 1.0
S  = SingleLayer(κ)
I  = Identity()
Dᵀ = DoubleLayerTransposed(κ)

@time Sxx  = assemble(S,X,X)
@time Dᵀxx = assemble(Dᵀ,X,X)
@time Ixx  = assemble(I,X,X)

E = PlaneWaveDirichlet(κ, point(1.0,0.0))
H = PlaneWaveNeumann(  κ, point(1.0,0.0))

Ex = assemble(E,X)
Hx = assemble(H,X)

x1 = Sxx \ Ex
x2 = (Ixx/2 + Dᵀxx) \ Hx

ℜ = real
ℑ = imag
n = numfunctions(X);
Δα = 2π/n;
α = (collect(1:n) - 0.5) * Δα
plt1 = plot(α, ℜ(x1), c=:blue, label="TM-EFIE")
scatter!(α, ℜ(x2), c=:red, m=:circle, label="TM-MFIE")
title!("current vs. angle")


# TE scattering
h = 2π / 51
Γ = meshcircle(1.0, h)

Y = lagrangec0d1(Γ)

N = HyperSingular(κ)
I = Identity()
D = DoubleLayer(κ)

@time Nyy = assemble(N,Y,Y)
@time Dyy = assemble(D,Y,Y)
@time Iyy = assemble(I,Y,Y)

E = PlaneWaveNeumann(  κ, point(1.0, 0.0))
H = PlaneWaveDirichlet(κ, point(1.0, 0.0))

Ey = assemble(E,Y)
Hy = assemble(H,Y)

x1 = Nyy \ Ey
x2 = (Iyy/2 - Dyy) \ Hy

n = numfunctions(Y);
Δα = 2π/n;
α = collect(0:n-1) * Δα

plt2 = plot(α, ℜ(x1), c=:blue, label="TE-EFIE")
scatter!(α, ℜ(x2), c=:red, m=:circle, label="TE-MFIE")
title!("current vs. angle")

plot(plt1,plt2)
