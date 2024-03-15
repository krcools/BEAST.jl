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

        if Mod.isstatic(gamma) # static case
            @assert !(isnothing(alpha)) && !(isnothing(beta))
        end

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

        if isnothing(alpha)
            if Mod.isstatic(gamma) # static case
                alpha = 1.0 # Default to double precision
            else
                alpha = one(gamma)
            end
        end

        Mod.MWDoubleLayer3D(alpha, gamma)
    end

    planewave(;
            direction    = error("missing arguement `direction`"),
            polarization = error("missing arguement `polarization`"),
            wavenumber   = error("missing arguement `wavenumber`"),
            amplitude    = one(real(typeof(wavenumber)))) =
        Mod.PlaneWaveMW(direction, polarization, wavenumber*im, amplitude)

    farfield(;
        wavenumber = error("missing argument: `wavenumber`")) =
            Mod.MWFarField3D(wavenumber=wavenumber)
end
