@testitem "GPU block assembly vs CPU" tags=[:gpu] begin
    using CUDA
    using CompScienceMeshes
    using Test

    @test CUDA.functional()
    CUDA.device!(parse(Int, get(ENV, "BEAST_GPU_DEVICE", "0")))

    extension = Base.get_extension(BEAST, :BEASTCUDAExt)
    mesh = meshcuboid(1.0, 1.0, 1.0, 0.5)
    sphere = readmesh(joinpath(@__DIR__, "..", "assets", "sphere5.in");
        T=Float64)
    quadrature = BEAST.DoubleNumSauterQstrat(2, 3, 2, 2, 2, 2)
    configurations = [
        (Helmholtz3D.singlelayer(wavenumber=1.0), lagrangec0(mesh; order=1)),
        (Maxwell3D.singlelayer(wavenumber=1.0), raviartthomas(mesh)),
        (Maxwell3D.singlelayer(wavenumber=1.0), buffachristiansen(sphere)),
    ]

    for (operator, space) in configurations
        cpu = assemble(operator, space, space;
            threading=:single, quadstrat=quadrature)
        assembler = extension.gpu_blockassembler(
            operator, space, space; quadstrat=quadrature)
        test_ids = collect(1:min(5, numfunctions(space)))
        trial_ids = collect(max(1, numfunctions(space) - 4):numfunctions(space))
        reverse!(trial_ids)
        workspace = CUDA.zeros(scalartype(operator, space, space),
            length(test_ids) + 2, length(trial_ids) + 2)
        block = @view workspace[2:end-1, 2:end-1]

        @test assembler(block, test_ids, trial_ids) === block
        @test isapprox(Array(block), cpu[test_ids, trial_ids];
            atol=1.0e-10, rtol=1.0e-10)
    end
end

@testitem "GPU far block assembly vs CPU" tags=[:gpu] begin
    using CUDA
    using CompScienceMeshes
    using Test

    @test CUDA.functional()
    CUDA.device!(parse(Int, get(ENV, "BEAST_GPU_DEVICE", "0")))

    extension = Base.get_extension(BEAST, :BEASTCUDAExt)
    test_mesh = meshcuboid(1.0, 1.0, 1.0, 0.5)
    trial_mesh = CompScienceMeshes.translate(
        test_mesh, point(5.0, 0.0, 0.0))
    test_space = lagrangec0(test_mesh; order=1)
    trial_space = lagrangec0(trial_mesh; order=1)
    operator = Helmholtz3D.singlelayer(wavenumber=1.0)
    quadrature = BEAST.DoubleNumQStrat(2, 3)
    test_ids = collect(1:4)
    trial_ids = collect(2:5)

    cpu = assemble(operator, test_space, trial_space;
        threading=:single, quadstrat=quadrature)
    assembler = extension.gpu_blockassembler(
        operator, test_space, trial_space; quadstrat=quadrature)
    block = CUDA.zeros(scalartype(operator, test_space, trial_space),
        length(test_ids), length(trial_ids))

    @test assembler(block, test_ids, trial_ids) === block
    @test isapprox(Array(block), cpu[test_ids, trial_ids];
        atol=1.0e-10, rtol=1.0e-10)
end
