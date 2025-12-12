@testitem "assemble linearmap" begin

    using CompScienceMeshes
    using LinearAlgebra

    fn = joinpath(pkgdir(BEAST), "test", "assets", "sphere45.in")
    m = CompScienceMeshes.readmesh(fn)

    X = raviartthomas(m)
    Id = BEAST.Identity()
    Ixx = BEAST.assemble(Id, X, X)
    IXX = BEAST.lu(Ixx)

    @test IXX isa BEAST.LinearMaps.LinearMap

    @hilbertspace m j
    @hilbertspace p q

    b = Ixx[p,m] + Ixx[q,j]
    a = IXX[p,m] + IXX[q,j]

    V = BEAST.DirectProductSpace([X,X])
    A = assemble(a, V, V)
    B = assemble(b, V, V)

    Ma = Matrix(A)
    Mb = Matrix(B)

    P = Matrix{Float64}(I, size(Ma,1), size(Mb,2))
    Q = Ma * Mb

    @test P ≈ Q atol=1e-10
end