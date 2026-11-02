@testitem "GWP interpolator" begin
    using CompScienceMeshes
    using LinearAlgebra

    dir = pkgdir(BEAST, "test/assets")
    Γ = CompScienceMeshes.readmesh(joinpath(dir, "sphere45.in"))
    ℰ = CompScienceMeshes.skeleton(Γ, 1)

    F = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=1.0)
    f = n × F

    X = BEAST.gwpdiv(Γ; order=3)
    Π = BEAST.GWPInterpolator(3, Γ, ℰ, Γ)

    u = Π(f)
    @test length(u) == numfunctions(X)

    # A quadratic field must be reproduced exactly by an order 2 interpolant
    struct QuadraticField end
    (::QuadraticField)(p) = (x = cartesian(p); point(x[1]*x[2], x[3]^2, x[1]-x[2]*x[3]))
    BEAST.scalartype(::QuadraticField) = Float64

    g = QuadraticField()
    Y = BEAST.gwpcurl(Γ, ℰ; order=2)
    v = BEAST.GWPInterpolator(2, Γ, ℰ, Γ)(g)
    ϕ = BEAST.refspace(Y)
    err = maximum(1:length(Γ)) do c
        p = neighborhood(chart(Γ, c), (0.2, 0.3))
        ϕp = ϕ(p)
        w = sum(v[i] * s.coeff * ϕp[s.refid].value
            for (i,fn) in enumerate(Y.fns) for s in fn if s.cellid == c)
        n̂ = normal(p)
        gp = g(p)
        norm(w - (gp - dot(n̂,gp)*n̂))
    end
    @test err < 1e-10
end