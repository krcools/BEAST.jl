using BEAST
using CompScienceMeshes
using LinearAlgebra
using Test
using StaticArrays

# define the plane wave excitation analytically
Γ = meshcuboid(1.0,1.0,1.0,0.2)
Y = unitfunctionc0d1(Γ)
X = raviartthomas(Γ)

### test crosstracemw
a = (@SVector [1.0,0.0,0.0])
eloc(x) = a*exp(-1.0im*dot((@SVector [1.0,0.0,0.0]),x))# not physical but to test the traces.
ex = BEAST.CrossTraceMW(BEAST.Trace(BEAST.FunctionWrapper{SVector{3,ComplexF64}}(eloc)))
# probeer testen te verzinnen die integraal nul geven.

rhs = assemble(ex, n× X)

# loop projector
∂Γ = boundary(Γ)
edges = setminus(skeleton(Γ,1), ∂Γ)
verts = setminus(skeleton(Γ,0), skeleton(∂Γ,0))
Λ = Matrix(connectivity(verts, edges, sign))

PΛ = Λ * pinv(Λ'*Λ) * Λ'

z = PΛ  * rhs
@test isapprox(norm(z), 0.0;atol = sqrt(eps()))

### test NDotTraceMW

a = (@SVector [1.0,0.0,0.0])
eloc(x) = a
ex = BEAST.NDotTrace(BEAST.Trace(BEAST.FunctionWrapper{SVector{3,Float64}}(eloc)))
# probeer testen te verzinnen die integraal nul geven.

rhs = assemble(ex, Y)
@test isapprox(norm(rhs), 0.0;atol = sqrt(eps()))


