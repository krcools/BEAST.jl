module Maxwell3D

    using ..BEAST
    Mod = BEAST

    """
        singlelayer(;gamma, alpha, beta)
        singlelayer(;wavenumber, alpha, beta)

    Bilinear form given by:

    ```math
        α ∬_{Γ×Γ} j(x)⋅k(y) G_{γ}(x,y) + β ∬_{Γ×Γ} div j(x) div k(y) G_{γ}(x,y)
    ```

    with ``G_{γ} = e^{-γ|x-y|} / 4π|x-y|``.
    """
    function singlelayer(;
            gamma=nothing,
            wavenumber=nothing,
            alpha=nothing,
            beta=nothing)

        
        gamma, wavenumber = Mod.gamma_wavenumber_handler(gamma, wavenumber)

        @assert gamma !== nothing

        alpha === nothing && (alpha = -gamma)
        beta  === nothing && (beta  = -1/gamma)

        Mod.MWSingleLayer3D(gamma, alpha, beta)
    end

    weaklysingular(;wavenumber) = singlelayer(;wavenumber, alpha=-im*wavenumber, beta=zero(im*wavenumber))
    hypersingular(;wavenumber) = singlelayer(; wavenumber, alpha=zero(im*wavenumber), beta=-1/(im*wavenumber))

    """
        doublelayer(;gamma)
        doublelaher(;wavenumber)

    Bilinear form given by:

    ```math
        α ∬_{Γ^2} k(x) ⋅ (∇G_γ(x-y) × j(y))
    ```

    with ``G_γ = e^{-γ|x-y|} / 4π|x-y|``
    """
    function doublelayer(;
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

        Mod.MWDoubleLayer3D(gamma)
    end
    """
        ndoublelayer(;gamma)
        ndoublelaher(;wavenumber)

    Bilinear form given by:

    ```math
        α ∬_{Γ^2} k(x) ⋅ (∇G_γ(x-y) × n_y j(y))
    ```

    with ``G_γ = e^{-γ|x-y|} / 4π|x-y|``
    """
    function ndoublelayer(;
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

        Mod.MWnDoubleLayer3D(gamma)
    end
    planewave(;
            direction    = error("missing arguement `direction`"),
            polarization = error("missing arguement `polarization`"),
            wavenumber   = error("missing arguement `wavenumber`"),
            amplitude    = one(real(typeof(wavenumber)))) =
        Mod.PlaneWaveMW(direction, polarization, wavenumber, amplitude)

    farfield(;
        wavenumber = error("missing argument: `wavenumber`")) =
            Mod.MWFarField3D(wavenumber=wavenumber)
    

    """
        ngreenint(;gamma)
        ngreenint(;wavenumber)

    Bilinear form given by:

    ```math
        α ∬_{Γ^2} k(x) ⋅ n_y (G_γ(x-y)j(y))
    ```

    with ``G_γ = e^{-γ|x-y|} / 4π|x-y|``
    """
    function ngreenint(;
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
        Mod.MWngreenint(gamma,alpha)
    end

    """
    greenint(;gamma)
    greenint(;wavenumber)

    Bilinear form given by:

    ```math
        α ∬_{Γ^2} k(x) ⋅ G_γ(x-y) j(y)
    ```

    with ``G_γ = e^{-γ|x-y|} / 4π|x-y|``
    """
    function greenint(;
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
        Mod.MWgreenint(gamma,alpha)
    end

    """
    greenint(;gamma)
    greenint(;wavenumber)

    Bilinear form given by:

    ```math
        α ∬_{Γ^2} k(x) ⋅ G_γ(x-y) j(y)
    ```

    with ``G_γ = e^{-γ|x-y|} / 4π|x-y|``
    """
    function greenint(;
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
        Mod.MWgreenint(gamma,alpha)
    end
    """
    gradgreenint(;gamma)
    gradgreenint(;wavenumber)

    Bilinear form given by:

    ```math
        α ∬_{Γ^2} k(x) ⋅∇G_γ(x-y) j(y)
    ```
    or
    ```math
    α ∬_{Γ^2} k(x) ∇G_γ(x-y) ⋅ j(y)
    ```
    with ``G_γ = e^{-γ|x-y|} / 4π|x-y|``
    """
    function gradgreenint(;
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
        Mod.MWgradgreenint(gamma,alpha)
    end
    """
    ngradgreenint(;gamma)
    ngradgreenint(;wavenumber)

    Bilinear form given by:

    ```math
        α ∬_{Γ^2} k(x) ⋅∇G_γ(x-y)n_y j(y)
    ```
    with ``G_γ = e^{-γ|x-y|} / 4π|x-y|``
    """
    function ngradgreenint(;
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
        Mod.MWngradgreenint(gamma,alpha)
    end
        
end
