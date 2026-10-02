@testitem "assemble BilForm: archive" begin

    using CompScienceMeshes
    using LinearAlgebra

    fn = joinpath(pkgdir(BEAST), "test", "assets", "sphere45.in")
    m = CompScienceMeshes.readmesh(fn)

    X = raviartthomas(m)
    𝕏 = X × X

    I = BEAST.Identity()
    T = Maxwell3D.singlelayer(wavenumber=1.0)

    @hilbertspace m j
    @hilbertspace k l

    a = (
        I[k,m] + 2*T[k,j] +
        3im*T[l,m] - I[l,j]
    )

    A = assemble(a, 𝕏, 𝕏; threading=:cellcoloring)
    # this test is brittle but not sure how to do this otherwise...
    @test A.maps[1].lmap.A.lmap === A.maps[4].lmap.A.lmap 
    @test A.maps[2].lmap.A.lmap === A.maps[3].lmap.A.lmap
end


@testitem "assemble BilForm: archive symmetric" begin

    using CompScienceMeshes
    using LinearAlgebra

    # fn = joinpath(pkgdir(BEAST), "test", "assets", "sphere45.in")
    # m = CompScienceMeshes.readmesh(fn)
    m = meshsphere(radius=1.0, h=0.6)

    RT = raviartthomas(m)
    BC = buffachristiansen(m)
    X = RT × BC

    # I = BEAST.Identity()
    T = Maxwell3D.singlelayer(wavenumber=1.0)
    K = Maxwell3D.doublelayer(wavenumber=1.0)

    @hilbertspace m j
    @hilbertspace k l

    a = (
        2*T[k,j] + 3im*T[l,m]
    )

    A = assemble(a, X, X; threading=:cellcoloring)
    @test A.maps[1].lmap.A.lmap === A.maps[2].lmap.A.lmap.parent 

    a = (
        2*K[k,j] + 3im*K[l,m]
    )

    A = assemble(a, X, X; threading=:cellcoloring)
    @test A.maps[1].lmap.A.lmap === A.maps[2].lmap.A.lmap.parent 
end