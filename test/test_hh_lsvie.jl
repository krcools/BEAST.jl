using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using SphericalScattering
using Test


# Homogeneous Dielectic Sphere Unit Test

@testset "Lippmann Schwinger Volume Integral Equation" begin

    # Environment
    ε1 = 1.0*SphericalScattering.ε0
    μ1 = SphericalScattering.μ0

    # Dielectic Sphere
    ε2 = 5.0*SphericalScattering.ε0
    μ2 = SphericalScattering.μ0
    r = 1.0
    sp = DielectricSphere(; radius = r, filling = Medium(ε2, μ2))


    # Mesh, Basis
    h = 0.35
    mesh = CompScienceMeshes.tetmeshsphere(r,h)
    bnd = boundary(mesh)
    X = lagrangec0d1(mesh; dirichlet = false)
    #@show numfunctions(X)
    strc = X -> strace(X,bnd)

    # VIE operators
    function generate_tau(ε_ins, ε_env)
        contr = 1.0 - ε_ins/ε_env
        function tau(x::SVector{U,T}) where {U,T}
            return T(contr)
        end
        return tau
    end
    τ = generate_tau(ε2, ε1)

    I, V, B =  Identity(), VIE.hhvolume(tau = τ, wavenumber = 0.0), VIE.hhboundary(tau = τ, wavenumber = 0.0)
    Y = VIE.hhvolumegradG(tau = τ, wavenumber = 0.0)


    # Exitation
    dirE = SVector(1.0, 0.0, 0.0)
    dirgradΦ = -dirE
    amp = 1.0

    # SphericalScattering exitation
    ex = UniformField(direction = dirE, amplitude = amp, embedding = Medium(ε1, μ1))

    # VIE exitation
    Φ_inc = VIE.linearpotential(direction=dirgradΦ, amplitude=amp)


    # Assembly
    b = real.(assemble(Φ_inc, X))

    Z_I = assemble(I, X, X)

    Z_V = assemble(V, X, X)
    Z_B = assemble(B, strc(X), X)

    Z_Y = assemble(Y, X, X)

    Z_version1 = Z_I + Z_V + Z_B
    Z_version2 = Z_I + Z_Y

    # MoM solution
    u_version1 = Z_version1 \ b
    u_version2 = Z_version2 \ b

    # Observation points
    range_ = range(-1.0*r,stop=1.0*r,length=14)
    points = [point(x,y,z) for x in range_ for y in range_ for z in range_]
    points_sp=[]
    for p in points
        norm(p)<0.97*r && push!(points_sp,p)
    end

    # SphericalScattering solution inside the dielectric sphere
    Φ = field(sp, ex, ScalarPotential(points_sp))

    # MoM solution inside the dielectric sphere
    Φ_MoM_version1 = BEAST.grideval(points_sp, u_version1, X)
    Φ_MoM_version2 = BEAST.grideval(points_sp, u_version2, X)



    err_Φ_version1 = norm(Φ - Φ_MoM_version1) / norm(Φ)
    err_Φ_version2 = norm(Φ - Φ_MoM_version2) / norm(Φ)

    @show err_Φ_version1
    @show err_Φ_version2

    @test err_Φ_version1 < 0.02
    @test err_Φ_version2 < 0.01

end