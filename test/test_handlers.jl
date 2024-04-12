using BEAST
using StaticArrays
using Test

for T in [Float16, Float32, Float64]
    @test_throws ErrorException BEAST.gamma_wavenumber_handler(T(1), T(2))

    gamma, wavenumber = BEAST.gamma_wavenumber_handler(T(1), nothing)
    @test gamma === T(1)

    gamma, wavenumber = BEAST.gamma_wavenumber_handler(nothing, T(2))
    @test gamma === T(2) * im

    gamma, wavenumber = BEAST.gamma_wavenumber_handler(nothing, nothing)
    @test BEAST.isstatic(gamma)
    @test gamma === Val(0)

    gamma, wavenumber = BEAST.gamma_wavenumber_handler(T(0), nothing)
    @test gamma === T(0)
    @test !BEAST.isstatic(gamma)

    gamma, wavenumber = BEAST.gamma_wavenumber_handler(nothing, T(1))
    @test gamma === T(1) * im
    @test !BEAST.isstatic(gamma)

    @test_throws ErrorException BEAST.operator_parameter_handler(T(1), T(1), T(1))
    alpha, gamma = BEAST.operator_parameter_handler(nothing, T(1), nothing)
    @test alpha === T(1)
    @test gamma === T(1)

    alpha, gamma = BEAST.operator_parameter_handler(nothing, nothing, nothing)
    @test alpha === 1.0
    @test BEAST.isstatic(gamma)

    alpha, gamma = BEAST.operator_parameter_handler(nothing, im * T(0), nothing)
    @test alpha === T(1)
    @test !BEAST.isstatic(gamma)

    # Helmholtz3D
    operator = Helmholtz3D.hypersingular(; alpha=nothing, beta=nothing, gamma=nothing,
        wavenumber=nothing
    )
    @test BEAST.isstatic(operator)
    @test operator.alpha === 0.0
    @test operator.beta === 1.0
    @test BEAST.gamma(operator) === 0.0

    operator = Helmholtz3D.hypersingular(; alpha=T(1), beta=T(2))
    @test BEAST.isstatic(operator)
    @test operator.alpha === T(1)
    @test operator.beta === T(2)
    @test BEAST.gamma(operator) === T(0)
    @test scalartype(operator) == T

    pwave = Helmholtz3D.planewave(; direction=SVector(T(0), T(0), T(1)))
    @test pwave.direction == SVector(T(0), T(0), T(1))
    @test pwave.gamma === T(0)
    @test pwave.amplitude === T(1)

    mpol = Helmholtz3D.monopole(; position=SVector(T(0), T(0), T(0)), amplitude=T(1))
    @test mpol.position === SVector(T(0), T(0), T(0))
    @test mpol.gamma === T(0)
    @test mpol.amplitude === T(1)

    gradmpol = Helmholtz3D.grad_monopole(; position=SVector(T(0), T(0), T(0)),
        amplitude=T(1))
    @test gradmpol.position === SVector(T(0), T(0), T(0))
    @test gradmpol.gamma === T(0)
    @test gradmpol.amplitude === T(1)

    # FarFields
    @test_throws AssertionError BEAST.MWFarField3D()
    @test_throws AssertionError BEAST.MWDoubleLayerFarField3D()
    @test_throws AssertionError BEAST.MWDoubleLayerRotatedFarField3D()

    # Maxwell3D
    operator = Maxwell3D.singlelayer(; alpha=T(1), beta=T(2))
    @test BEAST.isstatic(operator)
    @test operator.α === T(1)
    @test operator.β === T(2)
    @test BEAST.gamma(operator) === T(0)

    operator = Maxwell3D.doublelayer()
    @test BEAST.isstatic(operator)
    @test BEAST.gamma(operator) === 0.0

    operator = Maxwell3D.doublelayer(; alpha=T(1))
    @test BEAST.isstatic(operator)
    @test operator.alpha === T(1)
    @test BEAST.gamma(operator) === T(0)

    # BEAST
    pwave = BEAST.planewavemw3d(;
        direction=SVector(T(1), T(0), T(0)),
        polarization=SVector(T(0), T(1), T(0))
    )
    @test pwave.direction === SVector(T(1), T(0), T(0))
    @test pwave.polarisation === SVector(T(0), T(1), T(0))
    @test pwave.gamma === T(0)

    dpole = BEAST.dipolemw3d(location=SVector(T(0), T(0), T(0)),
        orientation=SVector(T(0), T(0), T(1)))
    @test dpole.location === SVector(T(0), T(0), T(0))
    @test dpole.orientation === SVector(T(0), T(0), T(1))
    @test dpole.gamma === T(0)
end
