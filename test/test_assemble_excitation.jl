@testitem "assemble exciation: multiple terms in linear form" begin

    using CompScienceMeshes
    using LinearAlgebra

    fn = joinpath(pkgdir(BEAST), "test", "assets", "sphere45.in")
    m = CompScienceMeshes.readmesh(fn)
    X = raviartthomas(m)

    κ = 1.0
    Einc1 = exp(+im*κ*0.5) * Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
    Einc2 = exp(-im*κ*0.5) * Maxwell3D.planewave(direction=-ẑ, polarization=-x̂, wavenumber=κ)
    
    e1 = (n × Einc1) × n
    e2 = (n × Einc2) × n

    @hilbertspace k
    l1 = e1[k]
    l2 = e2[k]
    l = e1[k] + e2[k]

    b = assemble(l, X)
    b1 = assemble(l1, X)
    b2 = assemble(l2, X)

    # @show norm(b)
    # @show norm(b1)
    # @show norm(b2)
    @test norm(b - (b1 + b2)) < 1e-10
end