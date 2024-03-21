module Helmholtz3D

    using ..BEAST
    const Mod = BEAST
    using LinearAlgebra

    function singlelayer(;
        alpha=nothing,
        gamma=nothing,
        wavenumber=nothing
    )

        alpha, gamma = Mod.operator_parameter_handler(alpha, gamma, wavenumber)

        return Mod.HH3DSingleLayerFDBIO(alpha,gamma)
    end

    function doublelayer(;
        alpha=nothing,
        gamma=nothing,
        wavenumber=nothing
    )

        alpha, gamma = Mod.operator_parameter_handler(alpha, gamma, wavenumber)

        return Mod.HH3DDoubleLayerFDBIO(alpha, gamma)
    end

    function doublelayer_transposed(;
        alpha=nothing,
        gamma=nothing,
        wavenumber=nothing
    )

        alpha, gamma = Mod.operator_parameter_handler(alpha, gamma, wavenumber)

        return Mod.HH3DDoubleLayerTransposedFDBIO(alpha, gamma)
    end

    function hypersingular(;
        alpha=nothing,
        beta=nothing,
        gamma=nothing,
        wavenumber=nothing
    )

        gamma, wavenumber = Mod.gamma_wavenumber_handler(gamma, wavenumber)

        if alpha === nothing
            if Mod.isstatic(gamma) #static case
                alpha = 0.0  # In the long run, this should probably be rather 'nothing'
            else
                alpha = gamma^2
            end

        end

        if beta === nothing
            if Mod.isstatic(gamma) #static case
                beta = one(alpha)
            else
                beta = one(gamma)
            end
        end

        return Mod.HH3DHyperSingularFDBIO(alpha, beta, gamma)
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

        return Mod.HH3DPlaneWave(direction, gamma, amplitude)
    end

    function linearpotential(; direction=SVector(1, 0, 0), amplitude=1.0)
        return Mod.HH3DLinearPotential(direction ./ norm(direction), amplitude)
    end

    function grad_linearpotential(; direction=SVector(0.0, 0.0, 0.0), amplitude=1.0)
        return Mod.gradHH3DLinearPotential(direction, amplitude)
    end

    function monopole(;
        position=SVector(0.0, 0.0, 0.0),
        gamma=nothing,
        wavenumber=nothing,
        amplitude=1.0
    )

        gamma, wavenumber = Mod.gamma_wavenumber_handler(gamma, wavenumber)
        Mod.isstatic(gamma) && (gamma = zero(amplitude))

        return Mod.HH3DMonopole(position, gamma, amplitude)
    end

    function grad_monopole(;
        position=SVector(0.0, 0.0, 0.0),
        gamma=nothing,
        wavenumber=nothing,
        amplitude=1.0
    )

        gamma, wavenumber = Mod.gamma_wavenumber_handler(gamma, wavenumber)
        Mod.isstatic(gamma) && (gamma = zero(amplitude))

        return Mod.gradHH3DMonopole(position, gamma, amplitude)
    end

end

export Helmholtz3D
