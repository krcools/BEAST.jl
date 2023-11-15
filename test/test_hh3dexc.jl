using CompScienceMeshes
using BEAST
using Test

for T in [Float32, Float64]
    sphere = readmesh(joinpath(dirname(@__FILE__),"assets","sphere5.in"), T=T)
    numcells(sphere)

    local κ = T(2π)
    direction = point(T,0,0,1)
    local f = BEAST.HH3DPlaneWave(direction, im*κ, T(1.0))

    v1 = f(point(T,0,0,0))
    v2 = f(point(T,0,0,0.5))

    @test v1 ≈ +1
    @test v2 ≈ -1

    lp = HH3DLinearPotential(point(T,0,1,0), 2.0)
    @test lp(point(T,1,1,0)) == T(2.0)

    gradlp = grad(lp)
    @test gradlp(point(T,1,1,0)) == point(T, 0, 2, 0)

    γnlp = dot(BEAST.NormalVector(), gradlp)

    import BEAST.∂n
    p = ∂n(f)

    local s = chart(sphere,first(sphere))
    local c = neighborhood(s, T.([1,1]/3))

    r = cartesian(c)
    local n = normal(s)

    w1 = p(c)
    w2 = -im*κ*dot(direction, n)*f(r)

    w1 ≈ w2

    local N = BEAST.HH3DHyperSingularFDBIO(im*κ)
    local X = BEAST.lagrangec0d1(sphere)

    numfunctions(X)

    BEAST.quadinfo(N, X, X)
    Nxx = assemble(N, X, X)

    @test size(Nxx) == (numfunctions(X), numfunctions(X))
end