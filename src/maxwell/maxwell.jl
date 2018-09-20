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

        if (gamma == nothing) && (wavenumber == nothing)
            error("Supply one of (not both) gamma or wavenumber")
        end

        if (gamma != nothing) && (wavenumber != nothing)
            error("Supply one of (not both) gamma or wavenumber")
        end

        if gamma == nothing
            if iszero(real(wavenumber))
                gamma = -imag(wavenumber)
            else
                gamma = im*wavenumber
            end
        end

        @assert gamma != nothing

        alpha == nothing && (alpha = -gamma)
        beta  == nothing && (beta  = -1/gamma)

        Mod.MWSingleLayer3D(gamma, alpha, beta)
    end

    """
        doublelayer(;gamma)
        doublelaher(;wavenumber)

    Bilinear form given by:

    ```math
        ∬_{Γ^2} k(x) ⋅ (∇G_γ(x-y) × j(y))
    ```

    with ``G_γ = e^{-γ|x-y|} / 4π|x-y|``
    """
    function doublelayer(;
            gamma=nothing,
            wavenumber=nothing)

        if (gamma == nothing) && (wavenumber == nothing)
            error("Supply one of (not both) gamma or wavenumber")
        end

        if (gamma != nothing) && (wavenumber != nothing)
            error("Supply one of (not both) gamma or wavenumber")
        end

        if gamma == nothing
            if iszero(real(wavenumber))
                gamma = -imag(wavenumber)
            else
                gamma = im*wavenumber
            end
        end

        @assert gamma != nothing

        Mod.MWDoubleLayer3D(gamma)
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
end
