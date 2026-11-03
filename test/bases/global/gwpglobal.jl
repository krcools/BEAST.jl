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
@testitem "GWP interpolator on spaces" begin
    using CompScienceMeshes
    using LinearAlgebra

    dir = pkgdir(BEAST, "test/assets")
    Γ = CompScienceMeshes.readmesh(joinpath(dir, "sphere45.in"))
    ℰ = CompScienceMeshes.skeleton(Γ, 1)

    # interpolating a GWP space on itself yields the identity
    Y2 = BEAST.gwpcurl(Γ, ℰ; order=2)
    Π2 = BEAST.GWPInterpolator(2, Γ, ℰ, Γ)
    A = Π2(Y2)
    @test size(A) == (numfunctions(Y2), numfunctions(Y2))
    @test norm(A - I, Inf) < 1e-10

    # lower order GWP functions are reproduced exactly
    Y1 = BEAST.gwpcurl(Γ, ℰ; order=1)
    B = Π2(Y1)
    @test size(B) == (numfunctions(Y2), numfunctions(Y1))

    function evaluate(X, coeffs, c, p)
        ϕp = BEAST.refspace(X)(p)
        sum(coeffs[i] * s.coeff * ϕp[s.refid].value
            for (i,fn) in enumerate(X.fns) for s in fn if s.cellid == c;
            init=zero(cartesian(p)))
    end

    err = maximum(1:numfunctions(Y1)) do m
        e_m = zeros(numfunctions(Y1)); e_m[m] = 1
        maximum(1:length(Γ)) do c
            p = neighborhood(chart(Γ, c), (0.2, 0.3))
            norm(evaluate(Y2, B[:,m], c, p) - evaluate(Y1, e_m, c, p))
        end
    end
    @test err < 1e-10
end

@testitem "GWP interpolator on spaces: consistency with embedding" begin
    using CompScienceMeshes
    using LinearAlgebra
    using SparseArrays

    dir = pkgdir(BEAST, "test/assets")
    Γ = CompScienceMeshes.readmesh(joinpath(dir, "sphere45.in"))
    ℰ = CompScienceMeshes.skeleton(Γ, 1)

    # a subset of the edges, some of which have their orientation flipped
    I = 1:2:length(ℰ)
    ℰ1 = Mesh(vertices(ℰ), [isodd(i ÷ 2) ? CompScienceMeshes.fliporientation(c) : c for (i,c) in zip(I, cells(ℰ)[I])])
    @test any(<(0), nonzeros(CompScienceMeshes.embedding(ℰ1, ℰ)))

    # for the lowest order, interpolation reduces to the signed injection of edges
    X1 = BEAST.gwpcurl(Γ, ℰ1; order=0)
    Π = BEAST.GWPInterpolator(0, Γ, ℰ, Γ)
    A = Π(X1)

    S = copy(transpose(CompScienceMeshes.embedding(ℰ1, ℰ)))
    @test size(A) == size(S)
    @test norm(A - S, Inf) < 1e-10
end
