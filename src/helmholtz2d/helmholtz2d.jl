module Helmholtz2D

    using ..BEAST
    const Mod = BEAST
    using LinearAlgebra

    function singlelayer(;
        alpha=nothing,
        gamma=nothing,
        wavenumber=nothing
    )

        alpha, gamma = Mod.operator_parameter_handler(alpha, gamma, wavenumber)

        return Mod.HH2DSingleLayerFDBIO(alpha, gamma)
    end

    function doublelayer(;
        alpha=nothing,
        gamma=nothing,
        wavenumber=nothing
    )

        alpha, gamma = Mod.operator_parameter_handler(alpha, gamma, wavenumber)

        return Mod.HH2DDoubleLayerFDBIO(alpha, gamma)
    end

    function doublelayer_transposed(;
        alpha=nothing,
        gamma=nothing,
        wavenumber=nothing
    )

        alpha, gamma = Mod.operator_parameter_handler(alpha, gamma, wavenumber)

        return Mod.HH2DDoubleLayerTransposedFDBIO(alpha, gamma)
    end

    function hypersingular(;
        alpha=nothing,
        gamma=nothing,
        wavenumber=nothing
    )
  
        alpha, gamma = Mod.operator_parameter_handler(alpha, gamma, wavenumber)

        return Mod.HH2DHyperSingularFDBIO(alpha, gamma)
    end

    function planewave(;
            direction=error("direction is a required argument"),
            gamma=nothing,
            wavenumber=nothing,
            amplitude=one(eltype(direction)))

        gamma, wavenumber = Mod.gamma_wavenumber_handler(gamma, wavenumber)

        # Note: Unlike for the operators, there seems little benefit in
        # explicitly declaring a Laplace-Type excitation.

        Mod.isstatic(gamma) && (gamma = zero(amplitude))
   
        return Mod.HH2DPlaneWave(direction, gamma, amplitude)
    end

end

export Helmholtz2D