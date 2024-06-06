using CompScienceMeshes
using BEAST

using Test

for U in [Float32, Float64]
    m1 = meshrectangle(U(1.0), U(1.0), U(0.5))
    m2 = CompScienceMeshes.translate!(meshrectangle(U(1.0), U(1.0), U(0.5)), point(U,0,0,1))
    m3 = CompScienceMeshes.translate!(meshrectangle(U(1.0), U(1.0), U(0.5)), point(U,0,0,2))

    X1 = raviartthomas(m1)
    X2 = raviartthomas(m2)
    X3 = raviartthomas(m3)

    local X = X1 × X2 × X3

    n1 = numfunctions(X1)
    n2 = numfunctions(X2)
    n3 = numfunctions(X3)

    @test numfunctions(X) == n1 + n2 + n3
    nt = numfunctions(X)

    T = MWSingleLayer3D(U(1.0))
    t = assemble(T, X, X)

    @test size(t) == (nt,nt)

    bilterms = [BEAST.Variational.BilTerm(1,1,Any[],Any[],1,T)]

    BilForm = BEAST.Variational.BilForm(:i, :j, bilterms)
    @test typeof(assemble(BilForm, X, X)) == LinearMaps.LinearCombination{U, Vector{LinearMap}}
end