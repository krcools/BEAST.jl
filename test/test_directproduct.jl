using CompScienceMeshes
using BEAST

using Test
import LinearMaps

for U in [Float32, Float64]
    local m1 = meshrectangle(U(1.0), U(1.0), U(0.5))
    local m2 = CompScienceMeshes.translate!(meshrectangle(U(1.0), U(1.0), U(0.5)), point(U,0,0,1))
    local m3 = CompScienceMeshes.translate!(meshrectangle(U(1.0), U(1.0), U(0.5)), point(U,0,0,2))

    local X1 = raviartthomas(m1)
    local X2 = raviartthomas(m2)
    local X3 = raviartthomas(m3)

    local X = X1 × X2 × X3

    local n1 = numfunctions(X1)
    local n2 = numfunctions(X2)
    local n3 = numfunctions(X3)

    @test numfunctions(X) == n1 + n2 + n3
    local nt = numfunctions(X)

    local T = MWSingleLayer3D(U(1.0))
    local t = assemble(T, X, X)

    @test size(t) == (nt,nt)

    local bilterms = [BEAST.Variational.BilTerm(1,1,Any[],Any[],1,T)]

    local BilForm = BEAST.Variational.BilForm(:i, :j, bilterms)
    local BXX = BEAST.assemble(BilForm, X, X)
    @assert BXX isa BEAST.LiftedMaps.LiftedMap
    # @test typeof(assemble(BilForm, X, X)) == LinearMaps.LinearCombination{U, Vector{LinearMaps.LinearMap}}
end