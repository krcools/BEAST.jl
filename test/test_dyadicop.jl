using Test

using LinearAlgebra
using CompScienceMeshes
using BEAST

@testset begin
    m = meshsphere(radius=1.0, h=0.35)
    X = lagrangec0d1(m)

    F = x -> 1.0
    f = BEAST.ScalarTrace{Float64}(F)
    g = BEAST.ScalarTrace{Float64}(F)

    op = 2*BEAST.DyadicOp(f,g)
    A = assemble(op, X, X)

    b = assemble(f,X)
    c = assemble(f,X)

    M = Matrix(A)
    @test norm(M - 2*b*c') ≤ 1e-8

    @hilbertspace j[1:2]
    @hilbertspace k[1:2]

    X2 = X × X
    A2 = assemble(@discretise(BEAST.diag(op)[k,j], k∈X2, j∈X2))
    M2 = Matrix(A2)

    nX = length(X)
    @test size(M2) == (2*nX, 2*nX)
end
