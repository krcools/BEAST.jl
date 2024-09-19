using CompScienceMeshes
using BEAST

using Test

@testset "Assembly of TensorOperator wrt SpaceTimeBasis" begin
    m = Mesh(
        [
            point(0,0,0),
            point(1,0,0),
            point(1,1,0),
            point(0,1,0)
        ],
        [
            index(1,2,3),
            index(1,3,4)
        ],
    )

    X = raviartthomas(m)
    δ = timebasisdelta(0.01234, 10)
    T = timebasisshiftedlagrange(0.1234, 10, 1)
    U = X ⊗ T
    V = X ⊗ δ

    Id = BEAST.Identity()
    # Nx = BEAST.NCross()
    op = Id ⊗ Id
    Z = assemble(op, V, U)
    A = BEAST.ConvolutionOperators.ConvOpAsArray(Z)
    @test size(A) == (1,1,2)
    for i in 2:10
        @test A[1,1,i] ≈ 0 atol=sqrt(eps(eltype(A)))
    end
end