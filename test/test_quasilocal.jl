using Test
using CompScienceMeshes, BEAST

@testitem "Quasi-local operator assembly" begin

    using CompScienceMeshes
    using LinearAlgebra
    using SparseArrays

    Γ = CompScienceMeshes.readmesh(joinpath(pkgdir(BEAST), "test", "assets", "sphere35.in"))
    RT = raviartthomas(Γ)
    Q = BEAST.QuasiLocalSingleLayerMW(1.0, 0.5)

    Z = assemble(Q, RT, RT; threading=:dofsplitting)
    @test Z isa SparseMatrixCSC
    @test 0 < nnz(Z) < length(Z) # the operator range truncates far interactions
    @test norm(Z) > 0

    # The storage is buffered per thread, so assemble from a task on every thread pool,
    # including the interactive one (thread numbering starts with the interactive pool).
    for pool in (:default, :interactive)
        Threads.nthreads(pool) == 0 && continue
        for threading in (:single, :dofsplitting)
            Zp = fetch(Threads.@spawn pool assemble(Q, RT, RT; threading))
            @test Zp ≈ Z
        end
    end
end
