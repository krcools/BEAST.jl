using Test
using LinearAlgebra
using CompScienceMeshes, BEAST

@testitem "integrate!: reusable quadrature buffer preserves fused integration" begin
    using BEAST
    using CompScienceMeshes
    using LinearAlgebra

    fn = joinpath(dirname(pathof(BEAST)), "../examples/assets/sphere45.in")
    m = BEAST.readmesh(fn)
    X = raviartthomas(m)
    op = Maxwell3D.singlelayer(gamma=1.0)

    # This quadrature strategy's `integrate!` method picks between five distinct
    # rule types (CommonFace/CommonEdge/CommonVertex/WiltonSERule/DoubleQuadRule)
    # for a given pair of triangles -- enough concrete types to defeat Julia's
    # union-splitting optimization if the chosen rule were allowed to escape as a
    # plain return value, instead of being consumed within the branch that built it.
    qs = BEAST.DoubleNumWiltonSauterQStrat(2, 3, 6, 7, 5, 5, 4, 3)

    tels, tad = assemblydata(X)
    bels, bad = assemblydata(X)
    trefs = brefs = refspace(X)
    qd = BEAST.quaddata(op, trefs, brefs, tels, bels, qs)
    qbuffer = BEAST.quadraturebuffer(qs, X, X)
    qaction = BEAST.ApplyIntegrate(qbuffer)

    zlocal = zeros(scalartype(op, X, X),
        numfunctions(trefs, CompScienceMeshes.domain(tels[1])),
        numfunctions(brefs, CompScienceMeshes.domain(bels[1])))

    buffered!(p, q) = begin
        tcell, bcell = tels[p], bels[q]
        fill!(zlocal, 0)
        BEAST.integrate!(op, trefs, brefs, p, tcell, q, bcell, qd, qs,
            zlocal, X, p, X, q; action=qaction)
    end

    # Exercise several branches rather than relying on a single quadrature rule.
    # This smoke check keeps the fused path bounded; the assembly test below is
    # the regression guard for worker-local buffer reuse.
    center(el) = sum(el.vertices) / length(el.vertices)
    p = 1
    distances = sortperm([norm(center(tels[p]) - center(bels[q])) for q in eachindex(bels)])
    test_qs = unique([distances[1], distances[2], distances[end]])

    for q in test_qs
        buffered!(p, q)
    end

    q = first(test_qs)
    @test_throws UndefKeywordError BEAST.integrate!(
        op,
        trefs,
        brefs,
        p,
        tels[p],
        q,
        bels[q],
        qd,
        qs,
        zlocal,
        X,
        p,
        X,
        q,
    )

    buffered_bytes = sum(q -> @allocated(buffered!(p, q)), test_qs)

    @test buffered_bytes < 10_000
    @test_throws MethodError BEAST.ApplyIntegrate()
    @test_throws MethodError BEAST.ApplyIntegrateNonConforming()
end

@testitem "assemblechunk_body! reuses Sauter-Schwab quadrature buffers" begin
    using BEAST
    using CompScienceMeshes
    using LinearAlgebra

    fn = joinpath(dirname(pathof(BEAST)), "../examples/assets/sphere45.in")
    m = BEAST.readmesh(fn)
    X = raviartthomas(m)
    op = Maxwell3D.singlelayer(gamma=1.0)
    qs = BEAST.DoubleNumWiltonSauterQStrat(2, 3, 6, 7, 5, 5, 4, 3)
    assembler = BEAST.blockassembler(op, X, X; quadstrat=qs)

    ids = collect(1:24)
    store(v, m, n) = nothing

    assembler(ids, ids, store)
    bytes = @allocated assembler(ids, ids, store)

    @test bytes < 3_000_000
end
