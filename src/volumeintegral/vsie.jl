module VSIE
    using ..BEAST
    Mod = BEAST

    """
        singlelayer(;gamma, alpha, beta, tau)
        singlelayer(;wavenumber, alpha, beta, tau)

    Bilinear form given by:

    ```math
        α ∬_{Γ×Ω} j(x)⋅τ(y)⋅k(y) G_{γ}(x,y) + β ∬_{Γ×Ω} div j(x) τ(y) div k(y) G_{γ}(x,y)
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

        alpha == nothing && (alpha = -gamma)
        beta  == nothing && (beta  = -1/gamma)
        tau   == nothing && (tau   = x->1.0)

        Mod.VSIESingleLayer(gamma, alpha, beta, tau)
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

        alpha == nothing && (alpha = -gamma)
        beta  == nothing && (beta  = -1/gamma)
        tau   == nothing && (tau   = x->1.0)

        Mod.VSIESingleLayer2(gamma, alpha, beta, tau)
    end

    """
    boundary(;gamma, alpha, tau)
    boundary(;wavenumber, alpha, tau)

    Bilinear form given by:

    ```math
        α ∬_{Γ×Ω} T_n j(x) grad G_{γ}(x,y)⋅τ(y) div k(y)
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

        Mod.VSIEBoundary(gamma, alpha, tau)
    end

    function boundaryT(;
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

        Mod.VSIEBoundaryT(gamma, alpha, tau)
    end

    

    """
    doublelayer(;gamma, alpha, beta, tau)
    doublelayer(;wavenumber, alpha, beta, tau)

    Bilinear form given by:

    ```math
        α ∬_{Ω×Γ} j(x) ⋅ grad G_{γ}(x,y) × τ(y)⋅k(y)
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

        Mod.VSIEDoubleLayer(gamma, alpha, tau)
    end


    function doublelayerT(;
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

        Mod.VSIEDoubleLayerT(gamma, alpha, tau)
    end

end