module Maxwell3D
    using ..BEAST
    Mod = BEAST



    """
    green(;gamma)
    green(;wavenumber)

    Bilinear form given by:

    ```math
        α ∬_{Γ^2} k(x) ⋅ G_γ(x-y) j(y)
    ```

    with ``G_γ = e^{-γ|x-y|} / 4π|x-y|``
    """
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
        Mod.HHHgreen(gamma,alpha)
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
    Mod.HHHgradgreen(gamma,alpha)
end
        
end