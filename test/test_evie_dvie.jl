using SphericalScattering
using BEAST, CompScienceMeshes
using StaticArrays
using LinearAlgebra


using Test


@testset "EVIE and DVIE - Scattering of plane waves on a sphere" begin 

    ϵ0 = 1.0
    ϵ1 = 5.0
    μ0 = 1.0
    ω = 1.0
    κ0 = ω * √(ϵ0*μ0)
    ϵ_r = ϵ1/ϵ0

    amplitude = 1.0
    polarization = x̂
    direction = ẑ

    E_inc = VIE.planewave(amplitude = amplitude, direction = direction, polarization = polarization, wavenumber = κ0)

    r = 1.0
    h = 0.35
    mesh = CompScienceMeshes.tetmeshsphere(r,h)
    bnd = boundary(mesh)
    #using PlotlyJS
    #PlotlyJS.plot(CompScienceMeshes.wireframe(mesh))


    # Analytical solution from SphericalScattering.jl
    sp = DielectricSphere(; radius = r, filling = Medium(ϵ1, μ0))
    ex = planeWave( amplitude = amplitude, 
                    polarization = polarization, 
                    direction = direction, 
                    embedding = Medium(ϵ0, μ0), 
                    frequency = ω/(2*pi))

    # Observation points
    range_ = range(-1.0*r,stop=1.0*r,length=24)
    points = [point(x,y,z) for x in range_ for y in range_ for z in range_]
    points_inside = Vector{SVector{3,Float64}}()
    for p in points
        norm(p)<0.95*r && push!(points_inside, p)
    end


    E_inside = field(sp, ex, ElectricField(points_inside))


    # evie
    function generate_tau2(ε_ins, ε_env)
        function tau(x::SVector{U,T}) where {U,T}
            return ε_ins/ε_env - 1.0        # (ϵ_r - 1)
        end
        return tau
    end
    tau2 = generate_tau2(ϵ1, ϵ0)

    ttrc = X2 -> ttrace(X2,bnd)

    X2 = nedelecc3d(mesh)
    @show numfunctions(X2)

    I2 = Identity()
    K2 = VIE.singlelayer2(wavenumber = κ0, tau = tau2)
    B2 = VIE.boundary2(wavenumber = κ0, tau = tau2)

    M_I2 = assemble(I2, X2, X2)
    M_K2 = assemble(-K2, X2, X2)
    M_B2 = assemble(B2, ttrc(X2), X2)

    M_evie = ϵ_r*M_I2 + M_K2 + M_B2

    b2 = assemble(E_inc, X2)
    u_evie = M_evie \ b2

    E_evie_inside = BEAST.grideval(points_inside, u_evie, X2)


    # dvie
    function generate_tau(ε_ins, ε_env)
        function tau(x::SVector{U,T}) where {U,T}
            return 1.0 - ε_env/ε_ins        # (1 - 1/ϵ_r)
        end
        return tau
    end
    tau = generate_tau(ϵ1, ϵ0)


    ntrc = X -> ntrace(X,bnd)

    X = nedelecd3d(mesh)
    @show numfunctions(X)

    I = Identity()
    K = VIE.singlelayer(wavenumber = κ0, tau = tau)
    B = VIE.boundary(wavenumber = κ0, tau = tau)

    M_I = assemble(I, X, X)
    M_K = assemble(-K, X, X) 
    M_B = assemble(-B, ntrc(X), X)

    M_dvie = (1.0/ϵ_r)*M_I + M_K/ϵ0 + M_B/ϵ0 # Greens function has the factor 1/(4πϵ₀), but ϵ₀ not in BEAST... => /ϵ0

    b = ϵ0 * assemble(E_inc, X)
    u_dvie = M_dvie \ b

    E_dvie_inside = (1/ϵ1)*(BEAST.grideval(points_inside, u_dvie, X)) # E_inside = (1/ϵ1) * D_inside



    
    err_evie = norm(E_evie_inside - E_inside) / norm(E_inside)
    err_dvie = norm(E_dvie_inside - E_inside) / norm(E_inside)
    err_ = norm(E_dvie_inside - E_evie_inside) / norm(E_evie_inside)

    @show err_evie
    @show err_dvie
    @show err_

    @test err_evie < 0.15
    @test err_dvie < 0.17
    @test err_ < 0.18
    
end