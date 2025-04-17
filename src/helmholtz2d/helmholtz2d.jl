module Helmholtz2D

    using ..BEAST
    const Mod = BEAST
    using LinearAlgebra

    function singlelayer(;
        #alpha=nothing,
        #gamma=nothing,
        wavenumber=nothing
    )

        #alpha, gamma = Mod.operator_parameter_handler(alpha, gamma, wavenumber)

        return Mod.SingleLayer(wavenumber)
    end

    function doublelayer(;
        #alpha=nothing,
        #gamma=nothing,
        wavenumber=nothing
    )

        #alpha, gamma = Mod.operator_parameter_handler(alpha, gamma, wavenumber)

        return Mod.DoubleLayer(wavenumber)
    end

    function doublelayer_transposed(;
        #alpha=nothing,
        #gamma=nothing,
        wavenumber=nothing
    )

        #alpha, gamma = Mod.operator_parameter_handler(alpha, gamma, wavenumber)

        return Mod.DoubleLayerTransposed(wavenumber)
    end

    function hypersingular(;
        #alpha=nothing,
        #beta=nothing,
        #gamma=nothing,
        wavenumber=nothing
    )
        #=
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
     =#
        return Mod.HyperSingular(wavenumber)
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