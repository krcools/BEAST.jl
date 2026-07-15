@testitem "hh2d mie higher order" begin
    using BEAST
    using CompScienceMeshes
    using StaticArrays
    using LinearAlgebra
    using Test
    using SpecialFunctions


    # We check only
    for order = [2; 3]

        ε0 = 8.854187821e-12
        μ0 = 4π*1e-7
        c0 = 1/sqrt(ε0*μ0)
        η0 = sqrt(μ0/ε0)

        f = 1e9 # 1 GHz
        ω = 2π*f
        λ = c0/f
        k = 2π/λ
        h = λ/10
        a = 1.0 # radius of the scatterer

        circle = CompScienceMeshes.meshcircle(a, h)
        X0 = lagrangecx(circle, order=order)
        X1 = lagrangec0(circle, order=order)

        # Computing the fields on a circle of radius rc
        rc = 500.0*a
        pts = meshcircle(rc, 0.3 * rc).vertices

        ## Analytical solutions

        dbesselj(n,x) = (besselj(n-1, x)  - besselj(n+1, x)) / 2
        dhankelh2(n,x) = (hankelh2(n-1, x)  - hankelh2(n+1, x)) / 2

        cart2polar(x,y) = SVector(sqrt(x^2 + y^2), atan(y, x))




    ## Analytical Solutions

        # TM planewave solution (based on Sec 6.4.2, Jin, Theory and Computation of Electromagnetic Fields])
        function TM_pec_planewave(E0, k, a, ρ, φ)
            n = 0

            # Jin (6.4.11)
            a_n(n) = -(1.0im)^(-n)*(besselj(n, k*a) / hankelh2(n, k*a))
            valEz(n) =  a_n(n) * hankelh2(n, k*ρ)*exp(im*n*φ)
            valHφ(n) = 1/(im * ω * μ0) * k * a_n(n) * dhankelh2(n,k*ρ) * exp(im*n*φ)
            valHρ(n) = -1/(im * ω * μ0)  * (im * n / ρ) * a_n(n) * hankelh2(n, k*ρ) * exp(im*n*φ)

            retEz = valEz(0)
            retHφ = valHφ(0)
            retHρ = valHρ(0)

            while true
                n += 1
                bufEz = valEz(n) + valEz(-n)
                if abs(bufEz) >= abs(retEz) * eps(eltype(real(retEz))) * 1e3
                    retEz += bufEz
                else
                    break
                end
            end

            n=0
            while true
                n += 1
                bufHφ = valHφ(n) + valHφ(-n)
                if abs(bufHφ) >= abs(retHφ) * eps(eltype(real(retHφ))) * 1e3
                    retHφ += bufHφ
                else
                    break
                end
            end

            n=0
            while true
                n += 1
                bufHρ = valHρ(n) + valHρ(-n)
                if abs(bufHρ) >= abs(retHρ) * eps(eltype(real(retHρ))) * 1e3
                    retHρ += bufHρ
                else
                    break
                end
            end

            return E0 * retEz, E0*(retHρ*SVector(cos(φ), sin(φ)) + retHφ*SVector(-sin(φ), cos(φ)))
        end

        function TM_pec_planewave_E(E0, k, a, pts)
            return [TM_pec_planewave(E0, k, a, cart2polar(p[1], p[2])...)[1] for p in pts]
        end

        function TM_pec_planewave_H(E0, k, a, pts)
            return [TM_pec_planewave(E0, k, a, cart2polar(p[1], p[2])...)[2] for p in pts]
        end


        # TE planewave solution
        function TE_pec_planewave(H0, k, a, ρ, φ)
            n = 0

            # Jin (6.4.19)
            b(n) = -(1.0im)^(-n)*dbesselj(n,k*a)/dhankelh2(n,k*a)

            # Hz is not needed, but maybe we use it later for MFIE
            valHz(n) = b(n) * hankelh2(n, k*ρ)*exp(im*n*φ)
            valEφ(n) = -1/(im * ω * ε0) * k * b(n) * dhankelh2(n, k*ρ)*exp(im*n*φ)
            valEρ(n) = 1/(im * ω * ε0)  * (im * n / ρ) * b(n) * hankelh2(n, k*ρ)*exp(im*n*φ)

            retHz = valHz(0)
            retEφ = valEφ(0)
            retEρ = valEρ(0)

            while true
                n += 1
                bufHz = valHz(n) + valHz(-n)
                if abs(bufHz) >= abs(retHz) * eps(eltype(real(retHz))) * 1e3
                    retHz += bufHz
                else
                    break
                end
            end

            n = 0
            while true
                n += 1
                bufEφ = valEφ(n) + valEφ(-n)
                if abs(bufEφ) >= abs(retEφ) * eps(eltype(real(retEφ))) * 1e3
                    retEφ += bufEφ
                else
                    break
                end
            end

            n = 0
            while true
                n += 1
                bufEρ = valEρ(n) + valEρ(-n)
                if abs(bufEρ) >= abs(retEρ) * eps(eltype(real(retEρ))) * 1e3
                    retEρ += bufEρ
                else
                    break
                end
            end

            return H0*retHz, H0*(retEρ*SVector(cos(φ), sin(φ)) + retEφ*SVector(-sin(φ), cos(φ)))
        end

        function TE_pec_planewave_H(H0, k, a, pts)
            return [TE_pec_planewave(H0, k, a, cart2polar(p[1],p[2])...)[1] for p in pts]
        end

        function TE_pec_planewave_E(H0, k, a, pts)
            return [TE_pec_planewave(H0, k, a, cart2polar(p[1],p[2])...)[2] for p in pts]
        end




    ## TM Scattering

        # TM-EFIE system matrix
        𝒮 = Helmholtz2D.singlelayer(; alpha=im*k*η0, wavenumber=k)
        M_TMEFIE = assemble(𝒮, X0, X0)

        # TM-MFIE system matrix
        𝒟ᵀ = Helmholtz2D.doublelayer_transposed(; wavenumber=k)
        Dᵀ = assemble(𝒟ᵀ, X0, X0)

        I0 = assemble(BEAST.Identity(), X0, X0)
        M_TMMFIE = +0.5*I0 + Dᵀ

        # 1. Excitation: Plane wave (The incoming planewave is along the x-axis)
        E0 = 1.0 # Amplitude
        Ez_pw_inc = Helmholtz2D.planewave(; amplitude=E0, wavenumber=k, direction=SVector(1.0, 0.0))

        ez_pw_inc = assemble(DirichletTrace(Ez_pw_inc), X0)
        j_TMEFIE_pw =  M_TMEFIE \ ez_pw_inc

        # TM-EFIE: We compute the scattered Ez component (scalar)
        Ez_pw_sca_num = -potential(HH2DSingleLayerNear(𝒮), pts, j_TMEFIE_pw, X0; type=ComplexF64)
        Ez_pw_sca_ana = TM_pec_planewave_E(E0, k, a, pts)

        # Will not improve until we have curvilinear elements
        # But must not become worse either
        # TODO: Once curvilinear elements are support: refine order of functions and
        # mesh in a lock step and update these values
        Ez_pw_sca_bounds = [0.002; 0.002; 0.0014]
        @test norm(Ez_pw_sca_ana - Ez_pw_sca_num) ./ norm(Ez_pw_sca_ana) < Ez_pw_sca_bounds[order]

        Ht_pw_sca_num = -potential(HH2DDoubleLayerTransposedNear(𝒟ᵀ), pts, j_TMEFIE_pw, X0; type=SVector{2,ComplexF64})
        Ht_pw_sca_ana = TM_pec_planewave_H(E0, k, a, pts)

        Ht_pw_sca_bounds = [0.002; 0.002; 0.0015]
        @test norm(Ref(ẑ) .× Ht_pw_sca_num - Ht_pw_sca_ana) ./ norm(Ht_pw_sca_ana) < Ht_pw_sca_bounds[order]

        # TM-MFIE
        Ht_pw_inc =  - 1.0 / (im * ω * μ0) * curl(Ez_pw_inc)
        ht_pw_inc = -assemble(TangentTrace(Ht_pw_inc), X0)

        j_TMMFIE_pw =  M_TMMFIE \ ht_pw_inc

        j_pw_sca_bounds = [0.05; 0.001; 0.0009]
        @test norm(j_TMMFIE_pw - j_TMEFIE_pw) / norm(j_TMEFIE_pw) < j_pw_sca_bounds[order]

            

    ## TE Scattering

        I = 1.0
        ρp = 2.0
        φp = 0.0
        r_inc = 1.0

        pts_inc =  meshcircle(r_inc, 0.6*r_inc).vertices

        # Scattered Field
        # TE-MFIE
        𝒟 = Helmholtz2D.doublelayer(; wavenumber=k)
        D = assemble(𝒟, X1, X1)
        I1 = assemble(BEAST.Identity(), X1, X1)
        M_TEMFIE = 0.5I1 - D

        # TE-EFIE
        𝒩 = Helmholtz2D.hypersingular(;alpha=im*ω*μ0, beta=1.0/(im*ω*ε0), wavenumber=k)
        M_TEEFIE = assemble(𝒩, X1, X1)

        # 1. Excitation: Planewave
        H0 = 1.0
        Hz_pw_inc = Helmholtz2D.planewave(; amplitude=H0, wavenumber=k, direction=SVector(1.0, 0.0))
        hz_pw_inc = assemble(DirichletTrace(Hz_pw_inc), X1)
        j_TEMFIE_pw = M_TEMFIE \ hz_pw_inc

        Hz_pw_sca_num = potential(HH2DDoubleLayerNear(𝒟), pts, j_TEMFIE_pw, X1; type=ComplexF64)
        Hz_pw_sca_ana = TE_pec_planewave_H(H0, k, a, pts)

        Hz_pw_sca_bounds = [0.0015; 0.0015; 0.0015]
        @test norm(Hz_pw_sca_num - Hz_pw_sca_ana) /norm(Hz_pw_sca_ana) < Hz_pw_sca_bounds[order]

        Et_pw_inc = 1 / (im * ω * ε0) * curl(Hz_pw_inc)
        @test Et_pw_inc.direction == Hz_pw_inc.direction
        @test Et_pw_inc.polarization == SVector(0.0, 1.0)
        @test Et_pw_inc.amplitude ≈ Et_pw_inc.gamma / (im * ω * ε0) atol=1e-15
        et_pw_inc = assemble(TangentTrace(Et_pw_inc), X1)

        j_TEEFIE_pw = M_TEEFIE \ et_pw_inc

        j_pw_sca_bounds = [0.007; 0.0002; 1.19e-5]
        @test norm(j_TEEFIE_pw - j_TEMFIE_pw)/norm(j_TEMFIE_pw) < j_pw_sca_bounds[order]

        Et_pw_sca_num = -potential(HH2DHyperSingularNear(𝒩), pts, j_TEEFIE_pw, X1; type=SVector{2, ComplexF64})
        Et_pw_sca_ana = TE_pec_planewave_E(H0, k, a, pts)

        # We compute the scattered Ez component (scalar)
        Et_pw_sca_bounds = [0.003; 0.0016; 0.0016]
        @test norm(Ref(ẑ) .× Et_pw_sca_num - Et_pw_sca_ana) / norm(Et_pw_sca_ana) < Et_pw_sca_bounds[order]

    end
end