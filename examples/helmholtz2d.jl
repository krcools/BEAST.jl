using CompScienceMeshes, BEAST

h = 2π / 51; Γ = meshcircle(1.0, h)
X, Y = lagrangecxd0(Γ), lagrangec0d1(Γ)

κ = 1.0
S = Helmholtz2D.singlelayer(; wavenumber=κ)
Dᵀ = Helmholtz2D.doublelayer_transposed(; wavenumber = κ)
D = Helmholtz2D.doublelayer(; wavenumber=κ)
N = Helmholtz2D.hypersingular(; wavenumber=κ)
I  = Identity()

# TM scattering
E = PlaneWaveDirichlet(κ, point(1.0,0.0))
H = PlaneWaveNeumann(  κ, point(1.0,0.0))

@hilbertspace u; @hilbertspace u′
@hilbertspace t; @hilbertspace t′

tm_efie = @discretise           S[u′,u] == E[u′] u∈X u′∈X
tm_mfie = @discretise (0.5I + Dᵀ)[t′,u] == H[t′] u∈X t′∈X

u1, u2 = solve(tm_efie), solve(tm_mfie)

# TE scattering
E = PlaneWaveNeumann(  κ, point(1.0, 0.0))
H = PlaneWaveDirichlet(κ, point(1.0, 0.0))

te_efie = @discretise          N[t,t′] == E[t′] t∈Y t′∈Y
te_mfie = @discretise (0.5I - D)[t,u′] == H[u′] t∈Y u′∈Y

t1, t2 = solve(te_efie), solve(te_mfie)


using Plots

nx = numfunctions(X);
Δα = 2π/nx; α = (collect(1:nx) .- 0.5) * Δα
plt1 = Plots.plot(α, real(u1), c=:blue, label="TM-EFIE")
Plots.scatter!(α, real(u2), c=:red, m=:circle, label="TM-MFIE")
Plots.title!("current vs. angle")

ny = numfunctions(Y);
Δα = 2π/ny; α = collect(0:ny-1) * Δα
plt2 = Plots.plot(α, real(t1), c=:blue, label="TE-EFIE")
Plots.scatter!(α, real(t2), c=:red, m=:circle, label="TE-MFIE")
Plots.title!("current vs. angle")

Plots.plot(plt1,plt2)
