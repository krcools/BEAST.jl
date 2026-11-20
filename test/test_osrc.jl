@testitem "OSRC mapping farfield" begin
    using BEAST, CompScienceMeshes
    using SphericalScattering
    using LinearAlgebra

    # start with the definiton of the needed physical properties:
    c = 2.99792458e8      #speed of light
    f = 2/2*c             #frequency
    κ = 2π * f / c        #wavenumber
    η = 1.0

    Γ = readmesh(joinpath(dirname(pathof(BEAST)), "..", "test", "assets", "sphere2.in"))
    X = raviartthomas(Γ)
    Nd = n×X

    # define now the operators and the excitation
    T = Maxwell3D.singlelayer(wavenumber=κ)
    Einc = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
    Hinc = -1/(im*κ*η) * curl(Einc)

    # map the trace of the scattered electric field onto the trace of the scattered magnetic field using the OSRC EtM map
    etm = BEAST.EtmOsrc(κ, 10, pi/2, 1.0)
    etm_map = assemble(etm, X, X)
    
    e = (n × Einc)
    exsc = -assemble(e, X)
    Gxx = assemble(BEAST.Identity(), X, X)
    iG = BEAST.lu(Gxx)
    exsc_primal = iG * exsc
    mapped_h_Nd_primal = etm_map * exsc_primal

    # add incident field to get total field
    h = (n × Hinc)
    h_Nd = assemble(h, Nd)
    h_Nd_primal = iG * h_Nd
    h_tot_OSRC = mapped_h_Nd_primal .+ h_Nd_primal

    # Compute the far field using the OSRC solution and compare with Mie series
    Θ, Φ = range(0.0,stop=π,length=24), 0.0
    points_cartFF = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]
    FF_OSRC = -im * f / (2 * c) * potential(MWFarField3D(T), points_cartFF, h_tot_OSRC, Nd)
    
    sp = PECSphere(radius=1.0)
    ex = planeWave(frequency=f)
    FF = scatteredfield(sp, ex, FarField(points_cartFF))

    error_tot = norm(FF_OSRC .- FF) / norm(FF)
    @test error_tot < 0.15
end

@testitem "Pade approximation" begin
    using BEAST

    # analytical coefficients check
    pade3 = BEAST.construct(BEAST.PadeApproximant(3, 0))
    @test pade3.Aj ≈ [0.053787171163038064, 0.1746458477080449, 0.27156698112891703]
    @test pade3.Bj ≈ [0.8117449009293668, 0.3887395330218429, 0.049515566048790455]
    @test pade3.R0 ≈ 7.0

    # test rotated
    Np, θ_p = 100, pi/3
    aj = [2/(2*Np+1)*sin(j*pi/(2*Np + 1))^2 for j in 1:Np]
    bj = [cos(j*pi/(2*Np + 1))^2 for j in 1:Np]
    RNp = z -> 1 + sum(aj[j]*z/(1+bj[j]*z) for j in 1:Np)
    RNp_rotated = z -> exp(-im*θ_p/2) * RNp((1+z) * exp(im*θ_p) - 1)
    z_values = [0, 1, 1000, im, -im, 1-im, 10-im, 100-im]
    f_rational = BEAST.scalarapproximation(BEAST.PadeApproximant(Np, θ_p))
    @test f_rational.(z_values) ≈ RNp_rotated.(z_values)

    # numerical convergency test in test points
    function pade_error(z, Np, θ_p)
        f_rational = BEAST.scalarapproximation(BEAST.PadeApproximant(Np, θ_p))
        return abs(f_rational(z) - sqrt(Complex(1+z)))
    end

    # test near but on the correct side of branch cut
    @test pade_error(1 - 0.1im, 10, pi/2) < 1e-3
    @test pade_error(1 - 0.1im, 100, pi/2) < 1e-10
    @test pade_error(0.0, 20, pi/2) < 1e-10
    @test pade_error(-10 - 0.1im, 20, pi/3) < 1e-3
    @test pade_error(-10 - 0.1im, 50, pi/3) < 1e-10
    @test pade_error(-1e4 - 0.1im, 500, pi/2) < 1e-3
    @test pade_error(-1e4 - 0.1im, 2000, pi/2) < 1e-10
end

@testitem "Inverse rational function" begin
    using BEAST

    approximant = BEAST.PadeApproximant(100, pi/2)
    fRational = BEAST.scalarapproximation(approximant)
    fInverseRational = BEAST.scalarapproximation(BEAST.InverseApproximant(approximant))

    z_values = [0, 1, 1000, im, -im, 1-im, 10-im, 100-im]
    value_fRational = 1 ./ fRational.(z_values)
    value_fInverseRational = fInverseRational.(z_values)

    @test value_fRational≈value_fInverseRational
end

@testitem "MtE and EtM check (higher order)" begin
    using BEAST, CompScienceMeshes
    using LinearAlgebra

    Γ = CompScienceMeshes.meshcuboid(1.0, 1.0, 1.0, 1.0)
    X =  BEAST.gwpdiv(Γ; order=1)
    G = assemble(BEAST.Identity(), X, X)

    κ = pi*1.0
    Np = 4
    mte = BEAST.MteOsrc(κ, Np, pi/2, 1.0; order=1)
    etm = BEAST.EtmOsrc(κ, Np, pi/2, 1.0; order=1)
    mte_map = assemble(mte, X, X)
    etm_map = assemble(etm, X, X)

    Idcheck1 = - G * etm_map * mte_map;
    Idcheck2 = - mte_map * G * etm_map;

    b = rand(numfunctions(X))
    x = Idcheck1*b
    @test norm(b - x) < 1e-8

    x = Idcheck2*b
    @test norm(b - x) < 1e-8
end