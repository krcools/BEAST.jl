module VIE
    using ..BEAST
    Mod = BEAST

    """
        singlelayer(;gamma, alpha, beta, tau)
        singlelayer(;wavenumber, alpha, beta, tau)

    Bilinear form given by:

    ```math
        α ∬_{Ω×Ω} j(x)⋅τ(y)⋅k(y) G_{γ}(x,y) + β ∬_{Ω×Ω} div j(x) τ(y)⋅k(y)⋅grad G_{γ}(x,y)
    ```

    with ``G_{γ} = e^{-γ|x-y|} / 4π|x-y|``
    and  ``τ(y)`` contrast dyadic
    """

    function singlelayer(;
            gamma=nothing,
            wavenumber=nothing,
            alpha=nothing,
            beta=nothing,
            tau=nothing)

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

        alpha == nothing && (alpha = wavenumber*wavenumber)
        beta  == nothing && (beta  = 1.0)
        tau   == nothing && (tau   = x->1.0)

        Mod.VIESingleLayer(gamma, alpha, beta, tau)
    end


    function singlelayer2(;
            gamma=nothing,
            wavenumber=nothing,
            alpha=nothing,
            beta=nothing,
            tau=nothing)

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

        alpha == nothing && (alpha = 1.0)
        beta  == nothing && (beta  = 1.0)
        tau   == nothing && (tau   = x->1.0)

        Mod.VIESingleLayer2(gamma, alpha, beta, tau)
    end

    """
    boundary(;gamma, alpha, tau)
    boundary(;wavenumber, alpha, tau)

    Bilinear form given by:

    ```math
        α ∬_{Ω×Ω} T_n j(x) grad G_{γ}(x,y)⋅τ(y)⋅k(y)
    ```

    with ``G_{γ} = e^{-γ|x-y|} / 4π|x-y|``
    and  ``τ(y)`` contrast dyadic
    """

    function boundary(;
            gamma=nothing,
            wavenumber=nothing,
            alpha=nothing,
            tau=nothing)

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

        alpha == nothing && (alpha = 1.0)
        tau   == nothing && (tau   = x->1.0)

        Mod.VIEBoundary(gamma, alpha, tau)
    end

    function boundary2(;
            gamma=nothing,
            wavenumber=nothing,
            alpha=nothing,
            tau=nothing)

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

        alpha == nothing && (alpha = 1.0)
        tau   == nothing && (tau   = x->1.0)

        Mod.VIEBoundary2(gamma, alpha, tau)
    end

    """
    doublelayer(;gamma, alpha, beta, tau)
    doublelayer(;wavenumber, alpha, beta, tau)

    Bilinear form given by:

    ```math
        α ∬_{Ω×Ω} j(x) ⋅ grad G_{γ}(x,y) × τ(y)⋅k(y)
    ```

    with ``G_{γ} = e^{-γ|x-y|} / 4π|x-y|``
    and  ``τ(y)`` contrast dyadic
    """

    function doublelayer(;
            gamma=nothing,
            wavenumber=nothing,
            alpha=nothing,
            beta=nothing,
            tau=nothing)

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

        alpha == nothing && (alpha = 1.0)
        tau   == nothing && (tau   = x->1.0)

        Mod.VIEDoubleLayer(gamma, alpha, tau)
    end

    planewave(;
        direction    = error("missing arguement `direction`"),
        polarization = error("missing arguement `polarization`"),
        wavenumber   = error("missing arguement `wavenumber`"),
        amplitude    = one(real(typeof(wavenumber)))) =
    Mod.PlaneWaveVIE(direction, polarization, wavenumber, amplitude)

    farfield(;
        wavenumber = error("missing argument: `wavenumber`"), 
        tau=error("missing argument: `tau`")) =
            Mod.VIEFarField3D(wavenumber=wavenumber, tau= tau)



    # Operators of the Lippmann Schwinger Volume Integral Equation,
    # which is based on the generalized Helmholtz Equation:
    """
        hhvolume(;gamma, alpha, tau)
        hhvolume(;wavenumber, alpha, tau)

    Bilinear form given by:

    ```math
        α ∬_{Ω×Ω} (grad j(x)) ⋅ G_{γ}(x,y)τ(y) (grad k(y))
    ```

    with ``G_{γ} = e^{-γ|x-y|} / 4π|x-y|``
    and  ``τ(y)`` contrast function
    """
    function hhvolume(;
            gamma=nothing,
            wavenumber=nothing,
            alpha=nothing,
            tau=nothing)

        if (gamma === nothing) && (wavenumber === nothing)
            error("Supply one of (not both) gamma or wavenumber")
        end

        if (gamma !== nothing) && (wavenumber !== nothing)
            error("Supply one of (not both) gamma or wavenumber")
        end

        if gamma === nothing
            if iszero(real(wavenumber))
                gamma = -imag(wavenumber)
            else
                gamma = im*wavenumber
            end
        end

        @assert gamma !== nothing

        alpha === nothing && (alpha = -1.0)
        tau   === nothing && (tau   = x->1.0)

        Mod.VIEhhVolume(gamma, alpha, tau)
    end


    """
    hhboundary(;gamma, alpha, tau)
    hhboundary(;wavenumber, alpha, tau)

    Bilinear form given by:

    ```math
        α ∬_{∂Ω×Ω} n̂(x) ⋅ j(x) G_{γ}(x,y) τ(y) (grad k(y))
    ```

    with ``G_{γ} = e^{-γ|x-y|} / 4π|x-y|``
    and  ``τ(y)`` contrast function
    and ``n̂(x)`` normal vector
    """
    function hhboundary(;
            gamma=nothing,
            wavenumber=nothing,
            alpha=nothing,
            tau=nothing)

        if (gamma === nothing) && (wavenumber === nothing)
            error("Supply one of (not both) gamma or wavenumber")
        end

        if (gamma !== nothing) && (wavenumber !== nothing)
            error("Supply one of (not both) gamma or wavenumber")
        end

        if gamma === nothing
            if iszero(real(wavenumber))
                gamma = -imag(wavenumber)
            else
                gamma = im*wavenumber
            end
        end

        @assert gamma !== nothing

        alpha === nothing && (alpha = 1.0)
        tau   === nothing && (tau   = x->1.0)

        Mod.VIEhhBoundary(gamma, alpha, tau)
    end


    """
    hhvolumek0(;gamma, alpha, tau)
    hhvolumek0(;wavenumber, alpha, tau)

    Bilinear form given by:

    ```math
        α ∬_{Ω×Ω} j(x) G_{γ}(x,y)τ(y) k(y)
    ```

    with ``G_{γ} = e^{-γ|x-y|} / 4π|x-y|``
    and  ``τ(y)`` contrast function
    """
    function hhvolumek0(;
        gamma=nothing,
        wavenumber=nothing,
        alpha=nothing,
        tau=nothing)

        if (gamma === nothing) && (wavenumber === nothing)
            error("Supply one of (not both) gamma or wavenumber")
        end

        if (gamma !== nothing) && (wavenumber !== nothing)
            error("Supply one of (not both) gamma or wavenumber")
        end

        if gamma === nothing
            if iszero(real(wavenumber))
                gamma = -imag(wavenumber)
            else
                gamma = im*wavenumber
            end
        end

        @assert gamma !== nothing

        alpha === nothing && (alpha = wavenumber*wavenumber)
        tau   === nothing && (tau   = x->1.0)

        Mod.VIEhhVolumek0(gamma, alpha, tau)
    end



    """
        hhvolumegradG(;gamma, alpha, tau)
        hhvolumegradG(;wavenumber, alpha, tau)

    Bilinear form given by:

    ```math
        α ∬_{Ω×Ω} j(x) grad_y(G_{γ}(x,y)) τ(y) ⋅ (grad k(y))
    ```

    with ``G_{γ} = e^{-γ|x-y|} / 4π|x-y|``
    and  ``τ(y)`` constant function
    """
    function hhvolumegradG(;
            gamma=nothing,
            wavenumber=nothing,
            alpha=nothing,
            tau=nothing)

        if (gamma === nothing) && (wavenumber === nothing)
            error("Supply one of (not both) gamma or wavenumber")
        end

        if (gamma !== nothing) && (wavenumber !== nothing)
            error("Supply one of (not both) gamma or wavenumber")
        end

        if gamma === nothing
            if iszero(real(wavenumber))
                gamma = -imag(wavenumber)
            else
                gamma = im*wavenumber
            end
        end

        @assert gamma !== nothing

        alpha === nothing && (alpha = -1.0)
        tau   === nothing && (tau   = x->1.0)

        Mod.VIEhhVolumegradG(gamma, alpha, tau)
    end

    linearpotential(;
    direction    = error("missing arguement `direction`"),
    amplitude    = one(real(typeof(direction[1])))) =
    Mod.LinearPotentialVIE(direction, amplitude)

end