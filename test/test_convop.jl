using Test
using BEAST
using CompScienceMeshes

Γ = meshrectangle(1.0, 1.0, 0.5, 3)
X = raviartthomas(Γ)
SL = derive(TDMaxwell3D.singlelayer(speedoflight=1.0))

Δt, Nt = 1.01, 10
δ = timebasisdelta(Δt,Nt)
h = timebasisc0d1(Δt,Nt)

SLh = (SL, X⊗δ, X⊗h)
A = assemble(SLh...)

@show size(A)