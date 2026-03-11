@testitem "assemble: <Nd, K Nd> - non-conforming" begin
    using CompScienceMeshes
    using LinearAlgebra

    m = readmesh(joinpath(pkgdir(BEAST), "test/assets/sphere45.in"))
    X = gradient(lagrangec0d1(m))
    Y = gradient(duallagrangec0d1(m))
    K = Maxwell3D.doublelayer(gamma=1.0)
    nxK = BEAST.DoubleLayerRotatedMW3D(1.0, K.gamma)

    # Note:
    # Nd = n x RT
    # RT = -n x Nd
    nxY = n × Y

    qs = BEAST.DoubleNumSauterQstrat{Int64, Int64}(2, 3, 4, 4, 4, 4)
    A = assemble(K, Y, X; quadstrat=qs)
    B = assemble(nxK, nxY, X; quadstrat=qs)

    @show norm(A)
    @show norm(B)

    @test norm(A-B) < 1e-10
    @test size(A) == (numfunctions(Y), numfunctions(X))
end