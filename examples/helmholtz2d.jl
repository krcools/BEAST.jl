using CompScienceMeshes, BEAST

h = 2π / 51; Γ = meshcircle(1.0, h)
X, Y = lagrangecxd0(Γ), lagrangec0d1(Γ)

κ = 1.0
S, Dᵀ = SingleLayer(κ), DoubleLayerTransposed(κ)
D, N  = DoubleLayer(κ), HyperSingular(κ)
I  = Identity()

# TM scattering
E = PlaneWaveDirichlet(κ, point(1.0,0.0))
H = PlaneWaveNeumann(  κ, point(1.0,0.0))

@hilbertspace u; @hilbertspace v

tm_efie = @discretise           S[v,u] == E[v] u∈X v∈X
tm_mfie = @discretise (0.5I + Dᵀ)[v,u] == H[v] u∈X v∈X

x1, x2 = solve(tm_efie), solve(tm_mfie)

# TE scattering
E = PlaneWaveNeumann(  κ, point(1.0, 0.0))
H = PlaneWaveDirichlet(κ, point(1.0, 0.0))

te_efie = @discretise          N[u,v] == E[v] u∈Y v∈Y
te_mfie = @discretise (0.5I - D)[u,v] == H[v] u∈Y v∈Y

y1, y2 = solve(te_efie), solve(te_mfie)


using Plots

nx = numfunctions(X);
Δα = 2π/nx; α = (collect(1:nx) - 0.5) * Δα
plt1 = plot(α, real(x1), c=:blue, label="TM-EFIE")
scatter!(α, real(x2), c=:red, m=:circle, label="TM-MFIE")
title!("current vs. angle")

ny = numfunctions(Y);
Δα = 2π/ny; α = collect(0:ny-1) * Δα
plt2 = plot(α, real(y1), c=:blue, label="TE-EFIE")
scatter!(α, real(y2), c=:red, m=:circle, label="TE-MFIE")
title!("current vs. angle")

plot(plt1,plt2)
