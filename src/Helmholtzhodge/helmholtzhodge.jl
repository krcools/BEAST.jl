module HHH
    using ..BEAST
    Mod = BEAST

    function green(;
            alpha=nothing,
            gamma=nothing,
            wavenumber=nothing)

        gamma, wavenumber = Mod.gamma_wavenumber_handler(gamma, wavenumber)

        if alpha === nothing
            if gamma !== nothing
                alpha = one(gamma)
            else
                alpha = 1.0 # Default to double precision
            end
        end
        Mod.HHHgreen(alpha,gamma,Mod.hhhidentity())
    end
 
    function gradgreen(;
        alpha=nothing,
        gamma=nothing,
        wavenumber=nothing)

    gamma, wavenumber = Mod.gamma_wavenumber_handler(gamma, wavenumber)

    if alpha === nothing
        if gamma !== nothing
            alpha = one(gamma)
        else
            alpha = 1.0 # Default to double precision
        end
    end
    Mod.HHHgradgreen(alpha,gamma,Mod.hhhidentity())
end
        
end