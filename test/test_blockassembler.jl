using BEAST
using CompScienceMeshes
using LinearAlgebra
using Test

@testset "blockassembler" begin
    r = 10.0
    lambda = 20 * r
    wavenumber = 2 * pi / lambda

    sphere = readmesh(joinpath(@__DIR__, "assets", "sphere5.in"); T=Float64)
    operator = Maxwell3D.doublelayer(; wavenumber)
    space = raviartthomas(sphere)

    matrix = assemble(operator, space, space)
    blockassembler = BEAST.blockassembler(operator, space, space)
    ids = collect(eachindex(space.fns))

    function assembleblock(assembler, rows=ids, columns=ids)
        out = zeros(ComplexF64, length(rows), length(columns))
        store(v, m, n) = (out[m, n] += v)
        assembler(rows, columns, store)
        return out
    end

    expected = assembleblock(blockassembler)
    @test matrix ≈ expected atol=eps(Float64)

    poolsize = length(blockassembler.scratchpool.available)
    scratchids = Set(objectid.(blockassembler.scratchpool.available))
    @test assembleblock(blockassembler) ≈ expected
    @test length(blockassembler.scratchpool.available) == poolsize
    @test Set(objectid.(blockassembler.scratchpool.available)) == scratchids

    failingstore(args...) = error("store failed")
    @test_throws ErrorException blockassembler(
        [first(ids)], [first(ids)], failingstore
    )
    @test length(blockassembler.scratchpool.available) == poolsize

    results = fetch.([
        Threads.@spawn assembleblock(blockassembler) for
        _ in 1:(2 * Threads.nthreads() + 1)
    ])
    @test all(result -> result ≈ expected, results)
    @test length(blockassembler.scratchpool.available) == poolsize
    @test Set(objectid.(blockassembler.scratchpool.available)) == scratchids
end
