using Test
using BEAST
using CompScienceMeshes

Γ = meshrectangle(1.0, 1.0, 0.5, 3)
X = raviartthomas(Γ)
fns = numfunctions(X)
id = Identity()⊗Identity()

Δt, Nt = 1.01, 10
δ = timebasisdelta(Δt,Nt)
h = timebasisc0d1(Δt,Nt)

idop = (id, X⊗δ, X⊗h)
G = assemble(idop...)
@test size(G) == (fns, fns, Nt)