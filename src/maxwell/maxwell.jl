module Maxwell3D

    using ..BEAST
    Mod = BEAST

    """
        singlelayer(;gamma, alpha, beta)
        singlelayer(;wavenumber, alpha, beta)

    Maxwell 3D single layer operator.

    Either gamma, or the wavenumber, or α and β must be provided.

    If α and β are not provided explitly, they are set to ``α = -γ`` and ``β = -1/γ`` with ``γ=\\mathrm{j} k``.
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

    """
        weaklysingular(;wavenumber)

    Weakly singular part of the Maxwell 3D single layer operator.

    ``α = -\\mathrm{j} k`` is set with the wavenumber ``k``.
    """
    weaklysingular(;wavenumber) = singlelayer(;wavenumber, alpha=-im*wavenumber, beta=zero(im*wavenumber))

    """
        hypersingular(;wavenumber)

    Hyper singular part of the Maxwell 3D single layer operator.

    ``β = -1/\\mathrm{j} k`` is set with the wavenumber ``k``.
    """
    hypersingular(;wavenumber) = singlelayer(; wavenumber, alpha=zero(im*wavenumber), beta=-1/(im*wavenumber))

    """
        doublelayer(;gamma)
        doublelayer(;wavenumber)

    Maxwell double layer operator.

    Either gamma or the wavenumber must be provided. Optionally, also alpha can be provided.

    If alpha is not provided explitly, it is set to ``α = 1``.
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

    """
        planewave(;
                direction    = error("missing arguement `direction`"),
                polarization = error("missing arguement `polarization`"),
                wavenumber   = error("missing arguement `wavenumber`"),
                amplitude    = one(real(typeof(wavenumber)))) 

    Time-harmonic plane wave.
    """
    planewave(;
            direction    = error("missing arguement `direction`"),
            polarization = error("missing arguement `polarization`"),
            wavenumber   = error("missing arguement `wavenumber`"),
            amplitude    = one(real(typeof(wavenumber)))) =
        Mod.PlaneWaveMW(direction, polarization, wavenumber*im, amplitude)


    planewaveExtractedKernel(;
            direction    = error("missing arguement `direction`"),
            polarization = error("missing arguement `polarization`"),
            wavenumber   = error("missing arguement `wavenumber`"),
            amplitude    = one(real(typeof(wavenumber)))) =
        Mod.PlaneWaveExtractedKernelMW(direction, polarization, wavenumber*im, amplitude)

    farfield(;
        wavenumber = error("missing argument: `wavenumber`")) =
            Mod.MWFarField3D(wavenumber=wavenumber)
end
