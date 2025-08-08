@testitem "storage: local operators" begin
    using CompScienceMeshes
    using SparseArrays
    for T in [Float32, Float64]
        fn = joinpath(@__DIR__, "assets/sphere5.in")
        m = readmesh(fn, T=T)

        Id = BEAST.Identity()
        X = BEAST.raviartthomas(m)

        Z1 = assemble(Id, X, X, storage_policy=Val{:densestorage})
        Z2 = assemble(Id, X, X, storage_policy=Val{:bandedstorage})
        Z3 = assemble(Id, X, X, storage_policy=Val{:sparsedicts})

        @test Z1 isa DenseMatrix
        @test Z2 isa SparseMatrixCSC
        @test Z3 isa SparseMatrixCSC

        @test Z1 ≈ Z2 atol=1e-8
        @test Z1 ≈ Z3 atol=1e-8
    end
end