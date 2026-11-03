@testitem "Lagrange interpolator" begin
    using CompScienceMeshes

    struct CubicField end
    (::CubicField)(p) = (x = cartesian(p); x[1]^3 - 2x[1]*x[2]^2 + x[3] + 1)
    BEAST.scalartype(::CubicField) = Float64

    f = CubicField()
    Γ = CompScienceMeshes.readmesh(joinpath(pkgdir(BEAST, "test/assets"), "sphere45.in"))
    nodes = skeleton(Γ, 0)
    edges = skeleton(Γ, 1)

    for order in 1:4
        X = BEAST.lagrangec0(Γ, nodes, edges; order)
        Π = BEAST.LagrangeInterpolator(order, Γ, nodes, edges, Γ)
        u = Π(f)
        @test length(u) == numfunctions(X)

        ϕ = BEAST.refspace(X)
        err = maximum(1:length(Γ)) do c
            p = neighborhood(chart(Γ, c), (0.2, 0.3))
            ϕp = ϕ(p)
            w = sum(u[i] * s.coeff * ϕp[s.refid].value
                for (i,fn) in enumerate(X.fns) for s in fn if s.cellid == c)
            abs(w - f(p))
        end

        # the cubic field is reproduced exactly from order 3 onwards
        order >= 3 && @test err < 1e-10
        order < 3 && @test err > 1e-6
    end
end
