@testitem "GPU batched block assembly vs CPU" tags=[:gpu] begin
    using CUDA
    using CompScienceMeshes
    using Test

    @test CUDA.functional()
    CUDA.device!(parse(Int, get(ENV, "BEAST_GPU_DEVICE", "0")))

    extension = Base.get_extension(BEAST, :BEASTCUDAExt)
    cuboid = meshcuboid(1.0, 1.0, 1.0, 0.5)
    sphere = readmesh(joinpath(@__DIR__, "..", "assets", "sphere5.in");
        T=Float64)
    quadrature = BEAST.DoubleNumSauterQstrat(2, 2, 2, 2, 2, 2)
    configurations = [
        (Helmholtz3D.singlelayer(wavenumber=1.0),
            lagrangec0(cuboid; order=1)),
        (Maxwell3D.singlelayer(wavenumber=1.0), raviartthomas(sphere)),
        (Maxwell3D.singlelayer(wavenumber=1.0), buffachristiansen(sphere)),
    ]

    for (operator, space) in configurations
        cpu = assemble(operator, space, space;
            threading=:single, quadstrat=quadrature)
        assembler = extension.gpu_blockassembler(
            operator, space, space; quadstrat=quadrature)
        num_ids = min(numfunctions(space), 24)
        first_ids = collect(1:div(num_ids, 2))
        last_ids = collect(div(num_ids, 2) + 1:num_ids)
        test_ids = [first_ids, first_ids, last_ids]
        trial_ids = [first_ids, last_ids, last_ids]

        blocks = extension.gpu_batched_blockassemble(
            assembler, test_ids, trial_ids)
        @test length(blocks) == length(test_ids)
        for block in eachindex(blocks)
            @test size(blocks[block]) ==
                (length(test_ids[block]), length(trial_ids[block]))
            @test isapprox(blocks[block],
                cpu[test_ids[block], trial_ids[block]];
                atol=1.0e-10, rtol=1.0e-10)
        end
    end
end

@testitem "GPU batched block assembly variants" tags=[:gpu] begin
    using CUDA
    using CompScienceMeshes
    using Test

    @test CUDA.functional()
    CUDA.device!(parse(Int, get(ENV, "BEAST_GPU_DEVICE", "0")))

    extension = Base.get_extension(BEAST, :BEASTCUDAExt)
    sphere = readmesh(joinpath(@__DIR__, "..", "assets", "sphere5.in");
        T=Float64)
    space = raviartthomas(sphere)
    operator = Maxwell3D.singlelayer(wavenumber=1.0)
    quadrature = BEAST.DoubleNumSauterQstrat(2, 2, 2, 2, 2, 2)
    cpu = assemble(operator, space, space;
        threading=:single, quadstrat=quadrature)
    assembler = extension.gpu_blockassembler(
        operator, space, space; quadstrat=quadrature)

    num_ids = numfunctions(space)
    test_ids = [Int[], [1], collect(1:min(num_ids, 20)),
        collect(1:min(num_ids, 3))]
    trial_ids = [Int[], [1], collect(1:min(num_ids, 20)),
        collect(4:min(num_ids, 6))]
    blocks = extension.gpu_batched_blockassemble(
        assembler, test_ids, trial_ids)

    for block in eachindex(blocks)
        @test size(blocks[block]) ==
            (length(test_ids[block]), length(trial_ids[block]))
        @test isapprox(blocks[block],
            cpu[test_ids[block], trial_ids[block]];
            atol=1.0e-10, rtol=1.0e-10)
    end

    test_ids = [collect(1:12), collect(13:24)]
    trial_ids = [collect(1:12), collect(13:24)]
    unchunked = extension.gpu_batched_blockassemble(
        assembler, test_ids, trial_ids)
    chunked = extension.gpu_batched_blockassemble(
        assembler, test_ids, trial_ids; budget=1)
    @test unchunked == chunked

    device_blocks = [CUDA.zeros(eltype(block), size(block))
        for block in unchunked]
    @test extension.gpu_batched_blockassemble!(
        device_blocks, assembler, test_ids, trial_ids) === device_blocks
    @test Array.(device_blocks) == unchunked

    assemblers = [assembler]
    @test extension.gpu_batched_blockassemble(
        assemblers, test_ids, trial_ids) == unchunked
end

@testitem "GPU batched far block assembly vs CPU" tags=[:gpu] begin
    using CUDA
    using CompScienceMeshes
    using Test

    @test CUDA.functional()
    CUDA.device!(parse(Int, get(ENV, "BEAST_GPU_DEVICE", "0")))

    extension = Base.get_extension(BEAST, :BEASTCUDAExt)
    test_mesh = meshcuboid(1.0, 1.0, 1.0, 0.5)
    trial_mesh = CompScienceMeshes.translate(
        test_mesh, point(10.0, 0.0, 0.0))
    test_space = lagrangec0(test_mesh; order=1)
    trial_space = lagrangec0(trial_mesh; order=1)
    operator = Helmholtz3D.singlelayer(wavenumber=1.0)
    quadrature = BEAST.DoubleNumQStrat(3, 4)
    cpu = assemble(operator, test_space, trial_space;
        threading=:single, quadstrat=quadrature)
    assembler = extension.gpu_blockassembler(
        operator, test_space, trial_space; quadstrat=quadrature)

    test_ids = [collect(1:10), collect(1:5)]
    trial_ids = [collect(1:10), collect(6:10)]
    blocks = extension.gpu_batched_blockassemble(
        assembler, test_ids, trial_ids)
    for block in eachindex(blocks)
        @test isapprox(blocks[block],
            cpu[test_ids[block], trial_ids[block]];
            atol=1.0e-10, rtol=1.0e-10)
    end
end

@testitem "GPU batched block assembly with multiple devices" tags=[:gpu] begin
    using CUDA
    using CompScienceMeshes
    using Test

    @test CUDA.functional()
    devices = collect(CUDA.devices())
    if length(devices) < 2
        @test_skip length(devices) >= 2
    else
        devices = devices[1:2]
        CUDA.device!(first(devices))

        extension = Base.get_extension(BEAST, :BEASTCUDAExt)
        mesh = meshcuboid(1.0, 1.0, 1.0, 0.5)
        space = lagrangec0(mesh; order=1)
        operator = Helmholtz3D.singlelayer(wavenumber=1.0)
        quadrature = BEAST.DoubleNumSauterQstrat(2, 2, 2, 2, 2, 2)
        cpu = assemble(operator, space, space;
            threading=:single, quadstrat=quadrature)
        assemblers = extension.gpu_blockassemblers(
            operator, space, space; quadstrat=quadrature, devices)

        num_ids = min(numfunctions(space), 16)
        middle = div(num_ids, 2)
        first_ids = collect(1:middle)
        last_ids = collect(middle + 1:num_ids)
        test_ids = [first_ids, first_ids, last_ids, last_ids]
        trial_ids = [first_ids, last_ids, first_ids, last_ids]
        blocks = extension.gpu_batched_blockassemble(
            assemblers, test_ids, trial_ids; budget=1)

        for block in eachindex(blocks)
            @test isapprox(blocks[block],
                cpu[test_ids[block], trial_ids[block]];
                atol=1.0e-10, rtol=1.0e-10)
        end
    end
end
