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
        (Maxwell3D.singlelayer(wavenumber=1.0), raviartthomas(mesh)),
    ]

    for (operator, space) in configurations
        cpu = assemble(operator, space, space;
            threading=:single, quadstrat=quadrature)
        gpu = assemble(operator, space, space;
            threading=:gpu, quadstrat=quadrature)

        @test isapprox(gpu, cpu; atol=1.0e-10, rtol=1.0e-10)
    end
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
