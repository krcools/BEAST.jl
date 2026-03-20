@info "Executing test_assemble_InverseOperator.jl"

using CompScienceMeshes
using LinearAlgebra
using BEAST
using Test

fn = joinpath(pkgdir(BEAST), "test", "assets", "sphere45.in")
Γ1 = CompScienceMeshes.readmesh(fn)
Γ2 = Mesh([point(x+3.0,y,z) for (x,y,z) in vertices(Γ1)], deepcopy(cells(Γ1)))
Γ = [Γ1, Γ2]
numdoms = 2

X = raviartthomas.(Γ)
Y = buffachristiansen.(Γ)

P = BEAST.DirectProductSpace([Xᵢ×Xᵢ for Xᵢ ∈ X])
Q = BEAST.DirectProductSpace([Yᵢ×Yᵢ for Yᵢ ∈ Y])

@hilbertspace m j
@hilbertspace k l

p = BEAST.hilbertspace(:p, numdoms)
q = BEAST.hilbertspace(:q, numdoms)

N = BEAST.NCross()

iNi = BEAST.InverseOperator(N, solver=BEAST.lu)[k, m] + BEAST.InverseOperator(N, solver=BEAST.GMRES)[l, j]
iNdiag = sum(iNi[q[i], p[i]] for i in 1:numdoms)

iA = BEAST.assemble(iNdiag, Q, P)

@test iA isa BEAST.LinearMaps.LinearMap

Ni = N[k, m] + N[l, j]
Ndiag = sum(Ni[q[i], p[i]] for i in 1:numdoms)

A = assemble(Ndiag, Q, P)
IA = BEAST.lu(Matrix(A))

@test IA isa BEAST.LinearMaps.LinearMap

@test Matrix(iA) ≈ Matrix(IA) atol=1e-6