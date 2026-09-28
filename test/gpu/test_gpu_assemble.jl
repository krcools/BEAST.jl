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
end
