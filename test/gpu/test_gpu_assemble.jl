@testitem "GPU dense assembly vs CPU" tags=[:gpu] begin
    using CUDA
    using CompScienceMeshes
    using Test

    @test CUDA.functional()
    CUDA.device!(parse(Int, get(ENV, "BEAST_GPU_DEVICE", "0")))

    mesh = meshcuboid(1.0, 1.0, 1.0, 0.5)
    quadrature = BEAST.DoubleNumSauterQstrat(2, 3, 2, 2, 2, 2)
    configurations = [
        (Helmholtz3D.singlelayer(wavenumber=1.0), lagrangec0(mesh; order=1)),
        (Helmholtz3D.singlelayer(wavenumber=1.0), duallagrangec0d1(mesh)),
        (Maxwell3D.singlelayer(wavenumber=1.0), raviartthomas(mesh)),
    ]

    for (operator, space) in configurations
        cpu = assemble(operator, space, space;
            threading=:single, quadstrat=quadrature)
        gpu = assemble(operator, space, space;
            threading=:gpu, quadstrat=quadrature)

        @test isapprox(gpu, cpu; atol=1.0e-10, rtol=1.0e-10)
    end

    empty_mesh = mesh[Int[]]
    product_space = BEAST.DirectProductSpace(
        [raviartthomas(mesh), raviartthomas(empty_mesh)])
    operator = Maxwell3D.singlelayer(wavenumber=1.0)
    cpu = assemble(operator, product_space, product_space;
        threading=:single, quadstrat=quadrature)
    gpu = assemble(operator, product_space, product_space;
        threading=:gpu, quadstrat=quadrature)
    @test Matrix(gpu) ≈ Matrix(cpu)

    trial_mesh = CompScienceMeshes.translate(mesh, point(5.0, 0.0, 0.0))
    test_space = lagrangec0(mesh; order=1)
    trial_space = lagrangec0(trial_mesh; order=1)
    operator = Helmholtz3D.singlelayer(wavenumber=1.0)
    regular_quadrature = BEAST.DoubleNumQStrat(2, 2)
    cpu = assemble(operator, test_space, trial_space;
        threading=:single, quadstrat=regular_quadrature)
    gpu = assemble(operator, test_space, trial_space;
        threading=:gpu, quadstrat=regular_quadrature)
    @test isapprox(gpu, cpu; atol=1.0e-10, rtol=1.0e-10)
end

@testitem "GPU mesh validation" tags=[:gpu] begin
    using CUDA
    using CompScienceMeshes
    using Test

    @test CUDA.functional()
    CUDA.device!(parse(Int, get(ENV, "BEAST_GPU_DEVICE", "0")))

    extension = Base.get_extension(BEAST, :BEASTCUDAExt)
    coarse_mesh = meshcuboid(1.0, 1.0, 1.0, 0.5)
    refined_mesh = barycentric_refinement(coarse_mesh)
    coarse_space = lagrangec0(coarse_mesh; order=1)
    refined_space = lagrangec0(refined_mesh; order=1)
    operator = Helmholtz3D.singlelayer(wavenumber=1.0)
    singular_quadrature = BEAST.DoubleNumSauterQstrat(2, 2, 2, 2, 2, 2)

    @test_throws ArgumentError assemble(operator, refined_space, coarse_space;
        threading=:gpu, quadstrat=singular_quadrature)
    @test_throws ArgumentError extension.gpu_blockassembler(
        operator, refined_space, coarse_space;
        quadstrat=singular_quadrature)

    planar_mesh = meshrectangle(1.0, 1.0, 0.5, 2)
    planar_space = lagrangec0d1(planar_mesh, skeleton(planar_mesh, 0))
    @test_throws ArgumentError assemble(operator, planar_space, planar_space;
        threading=:gpu, quadstrat=singular_quadrature)
    @test_throws ArgumentError extension.gpu_blockassembler(
        operator, planar_space, planar_space;
        quadstrat=singular_quadrature)
end

@testitem "GPU dense assembly with tiling" tags=[:gpu] begin
    using CUDA
    using CompScienceMeshes
    using Test

    @test CUDA.functional()
    CUDA.device!(parse(Int, get(ENV, "BEAST_GPU_DEVICE", "0")))

    extension = Base.get_extension(BEAST, :BEASTCUDAExt)
    mesh = meshcuboid(1.0, 1.0, 1.0, 0.5)
    operator = Helmholtz3D.singlelayer(wavenumber=1.0)
    space = lagrangec0(mesh; order=1)
    quadrature = BEAST.DoubleNumSauterQstrat(2, 3, 2, 2, 2, 2)
    tiling = extension.TilingStrategy(
        extension.EqualTiling(2), extension.EqualTiling(2))

    cpu = assemble(operator, space, space;
        threading=:single, quadstrat=quadrature)
    gpu = assemble(operator, space, space;
        threading=:gpu, quadstrat=quadrature,
        tilingstrat=tiling, nstreams=2)

    @test isapprox(gpu, cpu; atol=1.0e-10, rtol=1.0e-10)

    tiling = extension.TilingStrategy(
        extension.WorksizeTiling(5), extension.WorksizeTiling(7))
    gpu = assemble(operator, space, space;
        threading=:gpu, quadstrat=quadrature,
        tilingstrat=tiling, nstreams=3)
    @test isapprox(gpu, cpu; atol=1.0e-10, rtol=1.0e-10)

    @test_throws ArgumentError assemble(operator, space, space;
        threading=:gpu, quadstrat=quadrature, devices=Int[])
    @test_throws ArgumentError assemble(operator, space, space;
        threading=:gpu, quadstrat=quadrature,
        devices=[CUDA.device(), CUDA.device()])
end

@testitem "GPU dense assembly with multiple devices" tags=[:gpu] begin
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
        operator = Helmholtz3D.singlelayer(wavenumber=1.0)
        space = lagrangec0(mesh; order=1)
        quadrature = BEAST.DoubleNumSauterQstrat(2, 3, 2, 2, 2, 2)
        cpu = assemble(operator, space, space;
            threading=:single, quadstrat=quadrature)

        gpu = assemble(operator, space, space;
            threading=:gpu, quadstrat=quadrature,
            devices, nstreams=2)
        @test isapprox(gpu, cpu; atol=1.0e-10, rtol=1.0e-10)

        tiling = extension.TilingStrategy(
            extension.EqualTiling(2), extension.EqualTiling(2))
        gpu = assemble(operator, space, space;
            threading=:gpu, quadstrat=quadrature,
            devices, tilingstrat=tiling)
        @test isapprox(gpu, cpu; atol=1.0e-10, rtol=1.0e-10)
    end
end
