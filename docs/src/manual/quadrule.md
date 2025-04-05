# Designing your own quadrature rule and strategy

In the context of multi-trace solvers, testing and trial functions with logically separate but geometrically coinciding support can interact. In general these support can be equipped with completely indepenent meshes.

For meshes of flat faceted triangular panels, BEAST.jl defines the `BEAST.NonConformingIntegralOpQStrat` strategy. The constructor of `NonConformingIntegralOpQStrat` takes another quadratue strategy `ctrat` that is fit to deal with pairs of mutually conforming meshes of flat faceted meshes. For a pair comprising a test triangle and a trial triangle, the appropriate quadrature rule is chosen as follows:

- If the triangles are well-separated, `cstrat` is run on the pair and the resulting quadrature strategy is chosen.
- If the triangles have a single vertex in common, the Sauter-Schwab quadrature rule for common vertices is chosen.
- If the triangles overlap or have edges that overlap, they are both refined so that a geometrically conforming mesh is obtained. For each pair of triangles in this refinement, `cstrat` is run and the resulting quadrature strategy is chosen.

This quadrature strategy is powerful but requires a lot of geometric processing, resulting in significant increases in matrix assembly times.

When the user has additional knowledge about the precise geometric constellation of the interacting meshes, a more economic approach can be desirable. It is of course the user's responsibility not to use this economic approach in a context where it is not applicable.

In this tutorial we assume that the test mesh is conforming to the barycentric refinement of the trial mesh. We propose a quadrature strategy `NonConfTestBaryRefOfTrialQStrat`, parametrised by an underlying quadrature strategy fit for mutually conforming meshes, that implements the following algorithm:

- If the triangles are well-separated, `cstrat` is run on the pair and the resulting quadrature strategy is chosen.
- If the test and trial triangle share a vertex, the barycentric refinement of the trial triangle is constructed and the decision on the quadrature rule is deferred to `cstrat`.

Because of the a priori knowledge about the relative constellation of the meshes, the construction in the second case will automatically result in a pair of mutually conforming meshes that can safely be send off to `cstrat`.

The first step in defining a new quadrature strategy is the definition of the corresponding type:

```julia
struct NonConfTestBaryRefOfTrialQStrat{P} <: BEAST.AbstractQuadStrat
    conforming_qstrat::P
end
```

The semantics of the quadrature strategy are captured by the definition of a pair of methods for the functions **quaddata** and **quadrule**, respectively. The purpose of quaddata is the computation of cache data that can speed up the assembly. Typically, this involves computing and storing triangle normals, and the quadrature points and weights for the various quadrature rules that the quadrature strategy considers.

Because in the majority of cases, the computation of the interaction is deferred to the conforming quadrature strategy, the computation of the cache is forwarded to its method of quaddata, resulting simply in:

```julia
function BEAST.quaddata(a, X, Y, test_charts, trial_charts, 
    quadstrat::NonConfTestBaryRefOfTrialQStrat)

    return quaddata(a, X, Y, test_charts, trial_charts,
        quadstrat.conforming_qstrat)
end
```

The method of `quadrule` is key to the definition of the quadrature strategy and contains the the actual algorithm in charge of choosing the quadrature rule for any given pair of triangles:

```julia
function BEAST.quadrule(a, X, Y, i, test_chart, j, trial_chart, qd,
    quadstrat::NonConfTestBaryRefOfTrialQStrat)

    nh = BEAST._numhits(test_chart, trial_chart)
    nh > 0 && return TestInBaryRefOfTrialQRule(quadstrat.conforming_qstrat)
    return BEAST.quadrule(a, X, Y, i, test_chart, j, trial_chart, qd,
        quadstrat.conforming_qstrat)
end
```

The function body is essentially a one-to-one translation of the quadrature rule selection algorithm above to julia. The object `qd` passed to `quadrule` is the cache computed by quadrule. In our case, it is simply passed on to the underlying conforming quadrature rule in the case of well separated triangles.

Next, we define a type representing the quadrature rule that is used when test triangle and trial triangle are not well-separated:

```julia
struct TestInBaryRefOfTrialQRule{S}
    conforming_qstrat::S
end
```

The type `TestInBaryRefOfTrialQRule` refers to the quadrature rule that is responsible for the actual computation of the interactions in case of adjacency or overlap. Quadrature rules are implemented by specifying a method for the function `BEAST.momintegrals!`.

```julia
function BEAST.momintegrals!(out, op,
    test_functions, test_cell, test_chart,
    trial_functions, trial_cell, trial_chart,
    qr::TestInBaryRefOfTrialQRule)

    test_local_space = refspace(test_functions)
    trial_local_space = refspace(trial_functions)

    num_tshapes = numfunctions(test_local_space, domain(test_chart))
    num_bshapes = numfunctions(trial_local_space, domain(trial_chart))

    T = coordtype(test_chart)
    z, u, h, t = zero(T), one(T), T(1//2), T(1//3)

    c = CompScienceMeshes.point(T, t, t)
    v = (
        CompScienceMeshes.point(T, u, z),
        CompScienceMeshes.point(T, z, u),
        CompScienceMeshes.point(T, z, z))
    e = (
        CompScienceMeshes.point(T, z, h),
        CompScienceMeshes.point(T, h, z),
        CompScienceMeshes.point(T, h, h))

    X = (
        CompScienceMeshes.simplex(v[1], e[3], c),
        CompScienceMeshes.simplex(v[2], c, e[3]),
        CompScienceMeshes.simplex(v[2], e[1], c),
        CompScienceMeshes.simplex(v[3], c, e[1]),
        CompScienceMeshes.simplex(v[3], e[2], c),
        CompScienceMeshes.simplex(v[1], c, e[2]))

    C = CompScienceMeshes.cartesian(
        CompScienceMeshes.neighborhood(trial_chart, c))
    V = CompScienceMeshes.cartesian.((
        CompScienceMeshes.neighborhood(trial_chart, v[1]),
        CompScienceMeshes.neighborhood(trial_chart, v[2]),
        CompScienceMeshes.neighborhood(trial_chart, v[3])))
    E = CompScienceMeshes.cartesian.((
        CompScienceMeshes.neighborhood(trial_chart, e[1]),
        CompScienceMeshes.neighborhood(trial_chart, e[2]),
        CompScienceMeshes.neighborhood(trial_chart, e[3])))

    trial_charts = (
        CompScienceMeshes.simplex(V[1], E[3], C),
        CompScienceMeshes.simplex(V[2], C, E[3]),
        CompScienceMeshes.simplex(V[2], E[1], C),
        CompScienceMeshes.simplex(V[3], C, E[1]),
        CompScienceMeshes.simplex(V[3], E[2], C),
        CompScienceMeshes.simplex(V[1], C, E[2]))

    quadstrat = qr.conforming_qstrat
    qd = BEAST.quaddata(op, test_local_space, trial_local_space,
        (test_chart,), trial_charts, quadstrat)

    Q = zeros(T, num_tshapes, num_bshapes)
    out1 = zero(out)
    for (q,chart) in enumerate(trial_charts)
        qr1 = BEAST.quadrule(op, test_local_space, trial_local_space,
            1, test_chart, q ,chart, qd, quadstrat)
            
        BEAST.restrict!(Q, trial_local_space, trial_chart, chart, X[q])

        fill!(out1, 0)
        BEAST.momintegrals!(out1, op,
            test_functions, nothing, test_chart,
            trial_functions, nothing, chart, qr1)

        for j in 1:num_bshapes
            for i in 1:num_tshapes
                for k in 1:size(Q, 2)
                    out[i,j] += out1[i,k] * Q[j,k]
end end end end end
```

The algorithm constructs the charts of the barycentric refinement of the trial chart, which either share a vertex or an edge with the test chart, or completely coincide with the test chart. Because of this, contributions from any of the refinement charts can be computed accuratey by a classic quadrature rule:

```julia
momintegrals!(outq, op,
    test_functions, nothing, test_chart,
    trial_functions, nothing, chart, qr)
```

This call to `momintegrals!` calculates interactions between the shape functions on the original test chart and the shape functions on one of the six subcharts in the refinement of the trial chart. The corresponding contribution to the interaction with the shape functions on the coarse trial chart can be calculated if we know how the restriction of the coarse shape functions to any of the subcharts can be written as linear combinations of the shape on that subchart. This information is provided by

```julia
BEAST.restrict!(Q, trial_local_space, trial_chart, chart, X[q])
```

For efficiency, the overlap function from the domain of `chart` to the domain of `trial_chart` has to be supplied.

!!! note
    The quadrature strategy and related quadrature rules implemented here can be rearded as meta-strategies, and meta-rules, as they defer most of the heavy lifting to underlying strategies and rules for mutually conforming meshes.

    In a *primitive* rule, methods of `momintegrals!` typically contain implementations of numerical quadrature methods.

To verify correctness of the above strategy, we can compare the results against existing routines that either provide less accurate results or similar results at reduced efficiency.

```julia
@testitem "NonConfTestBaryRefOfTrialQStrat" begin
    using CompScienceMeshes

    fnm = joinpath(dirname(pathof(BEAST)), "../test/assets/sphere45.in")
    Γ1 = BEAST.readmesh(fnm)
    Γ2 = deepcopy(Γ1)

    X = raviartthomas(Γ1)
    Y1 = buffachristiansen(Γ1)
    Y = buffachristiansen(Γ2)

    K = Maxwell3D.doublelayer(gamma=1.0)
    qs1 = BEAST.DoubleNumWiltonSauterQStrat(2, 3, 6, 7, 5, 5, 4, 3)
    qs2 = BEAST.NonConformingIntegralOpQStrat(qs1)
    qs3 = BEAST.NonConfTestBaryRefOfTrialQStrat(qs1)

    @time Kyx1 = assemble(K, Y, X; quadstrat=qs1)
    @time Kyx2 = assemble(K, Y, X; quadstrat=qs2)
    @time Kyx3 = assemble(K, Y, X; quadstrat=qs3)
    @time Kyx4 = assemble(K, Y1, X; quadstrat=qs1)

    using LinearAlgebra
    @test norm(Kyx1 - Kyx2) < 0.05
    @test norm(Kyx1 - Kyx3) < 0.05
    @test norm(Kyx2 - Kyx3) < 0.002
    @test norm(Kyx2 - Kyx4) < 0.002
    @test norm(Kyx3 - Kyx4) < 1e-12
end
```
