using CompScienceMeshes
using BEAST
using LinearAlgebra

width, height, h = 1.0, 0.5, 0.1
ℑ₁ = meshrectangle(width, height, h)
ℑ₂ = CompScienceMeshes.rotate(ℑ₁, 0.5π * x̂)
ℑ₃ = CompScienceMeshes.rotate(ℑ₁, 1.0π * x̂)
Γ₁ = weld(ℑ₁,-ℑ₂)
Γ₂ = weld(ℑ₂,-ℑ₃)

X₁ = raviartthomas(Γ₁)
X₂ = raviartthomas(Γ₂)

Y₁ = buffachristiansen(Γ₁)
Y₂ = buffachristiansen(Γ₂)

X = X₁ × X₂
Y = Y₁ × Y₂

@hilbertspace p
@hilbertspace q

κ = 3.0
SL = Maxwell3D.singlelayer(wavenumber=κ)
N = NCross()
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n;

import BEAST: diag, blocks
Sxx = assemble(@discretise blocks(SL)[p,q] p∈X q∈X)
ex = assemble(e, X)

# Build the preconditioner
Nxy = assemble(@discretise(BEAST.diag(N)[p,q], p∈X, q∈Y))
Dyx = BEAST.GMRESSolver(Nxy, tol=2e-12, restart=250, verbose=false)
Dxy = BEAST.GMRESSolver(transpose(Nxy), tol=2e-12, restart=250, verbose=false)
Syy = BEAST.assemble(@discretise BEAST.diag(SL)[p,q] p∈Y q∈Y)
P = Dxy * Syy * Dyx

# Solve system without and with preconditioner
u1, ch1 = solve(BEAST.GMRESSolver(Sxx,tol=2e-5, restart=250), ex)
u2, ch2 = solve(BEAST.GMRESSolver(P*Sxx, tol=2e-5, restart=250), P*ex)

# Compute and visualise the far field
Φ, Θ = [0.0], range(0,stop=π,length=100)
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]
near1 = potential(MWFarField3D(wavenumber=κ), pts, u1, X)
near2 = potential(MWFarField3D(wavenumber=κ), pts, u2, X)

using Plots
@show norm(u1-u2)
plot(title="far field")
plot!(Θ, real.(getindex.(near1,1)), label="no preconditioner")
scatter!(Θ, real.(getindex.(near2,1)), label="Calderon preconditioner")