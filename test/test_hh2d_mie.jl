using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test
using SpecialFunctions

const ε0 = 8.854187821e-12
const μ0 = 4π*1e-7
const c0 = 1/sqrt(ε0*μ0)
const η0 = sqrt(μ0/ε0)

f = 1e9 # 1 GHz
ω = 2π*f
λ = c0/f
k = 2π/λ
h = λ/10
a = 1.0 # radius of the scatterer

circle = CompScienceMeshes.meshcircle(a, h)
X0 = lagrangecxd0(circle)
X1 = lagrangec0d1(circle)

# Computing the fields on a circle of radius rc
rc = 500.0*a
pts = meshcircle(rc, 0.3 * rc).vertices

## Analytical solutions

dbesselj(n,x) = (besselj(n-1, x)  - besselj(n+1, x)) / 2
dhankelh2(n,x) = (hankelh2(n-1, x)  - hankelh2(n+1, x)) / 2

cart2polar(x,y) = SVector(sqrt(x^2 + y^2), atan(y, x))

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

testez(x, y) = TM_pec_planewave_E(1.0, k, a, [SVector(x, y)])[1]
testht(x, y) = TM_pec_planewave_H(1.0, k, a, [SVector(x, y)])[1]

function fdcurl(F::Function, p, h)
    x = p[1]
    y = p[2]
    curlF_x = (F(x, y+h) - F(x, y-h))/(2h)
    curlF_y = -(F(x+h, y) - F(x-h, y))/(2h)

    return SVector(curlF_x, curlF_y)
end

p = SVector(3.0, 2.0)

@test norm(-1/(im * ω * μ0) .* fdcurl(testez, p, 0.00001) - testht(p...)) ./ norm(testht(p...)) < 1e-8

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

testhz(x, y) = TE_pec_planewave_H(1.0, k, a, [SVector(x, y)])[1]
testet(x, y) = TE_pec_planewave_E(1.0, k, a, [SVector(x, y)])[1]

@test norm(1/(im * ω * ε0) .* fdcurl(testhz, p, 0.00001) - testet(p...)) ./ norm(testet(p...)) < 1e-8

"""
 Here, the line current is infinitely long along the z-axis, and is pointed towards the 
 z-axis as well. Thus, the potentials Ax, Ay will not be excited, leaving us only Az to 
 compute, which will then give us Ez. The monopole excitation in BEAST is used for this 
 excitation.

Ref: Sec 6.5.3, Jin, Theory and Computation of Electromagnetic fields
Analytical expression used for scattered field: Eq 6.5.38
"""
function TM_pec_line_curr(I, k, a, ρ, φ, ρp, φp)
    n = 0
    A_d(n) = besselj(n,k*a)/hankelh2(n,k*a)
    val(n) = A_d(n) * hankelh2(n,k*ρ) * hankelh2(n,k*ρp) * exp(im*n*(φ-φp))
    ret = val(0)
    while true
        n += 1
        buf = val(n) + val(-n)
        if abs(buf) >= abs(ret) * eps(eltype(real(ret))) * 1e3
            ret += val(n) + val(-n)
        else
            break
        end
    end
    return η0*k*I/4 * ret
end

function TM_pec_line_curr_E(I, k, a, pts, spt)
    cart2polar(x,y) = SVector(sqrt(x^2 + y^2), atan(y, x))
    return [TM_pec_line_curr(I,k,a,cart2polar(p[1],p[2])..., spt[1],spt[2]) for p in pts]
end


## TM scattering

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
Ez_pw_sca_num = -im*k*η0 * potential(HH2DSingleLayerNear(im * k), pts, j_TMEFIE_pw, X0; type=ComplexF64)
Ez_pw_sca_ana = TM_pec_planewave_E(E0, k, a, pts)

norm(Ez_pw_sca_ana - Ez_pw_sca_num) ./ norm(Ez_pw_sca_ana)

Ht_pw_sca_num = -potential(HH2DDoubleLayerTransposedNear(im * k), pts, j_TMEFIE_pw, X0; type=SVector{2,ComplexF64})
Ht_pw_sca_ana = TM_pec_planewave_H(E0, k, a, pts)

norm(Ht_pw_sca_num - Ht_pw_sca_ana) ./ norm(Ht_pw_sca_ana)

-1/(im * ω * μ0) .* fdcurl(testez, pts[end], 0.00001)
Ht_pw_sca_num[end]
Ht_pw_sca_num[4]
Ht_pw_sca_ana[4]
# TM-MFIE
Ht_pw_inc =  - 1.0 / (im * ω * μ0) * curl(Ez_pw_inc)
ht_pw_inc = assemble(TangentTrace(Ht_pw_inc), X0)

j_TMMFIE_pw =  M_TMMFIE \ ht_pw_inc

@test norm(j_TMMFIE_pw - j_TMEFIE_pw) / norm(j_TMEFIE_pw) < 0.05

# 2. Excitation: Infinite line current

I = 1.0 # 1 A
ρp = 2.0 # Source radial pos
φp = 0.0 # Source angular pos
xp = ρp * cos(φp)
yp = ρp * sin(φp)

# Choosing the amplitude of Einc such that [Eq 6.5.11, Jin] is satisfied
Ez_lc_Einc = Helmholtz2D.monopole(;position = SVector(xp,yp), wavenumber=k, amplitude = -η0*k*I/4)
ez_lc_Einc = assemble(DirichletTrace(Ez_lc_Einc),X0)
j_TMEFIE_lc = M_TMEFIE \ ez_lc_Einc

Ez_lc_sca_num = -im*η0*k * potential(HH2DSingleLayerNear(im*k),pts , j_TMEFIE_lc, X0;type=ComplexF64)

Ez_lc_sca_ana = TM_pec_line_curr_E(I, k, a, pts, SVector(ρp,φp))

@test norm(Ez_lc_sca_num - Ez_lc_sca_ana) / norm(Ez_lc_sca_ana) < 0.003


## TE scattering

# TE-MFIE
𝒟 = Helmholtz2D.doublelayer(; wavenumber=k)
D = assemble(𝒟, X1, X1)
I1 = assemble(BEAST.Identity(), X1, X1)
M_TEMFIE = 0.5I1 - D

# TE-EFIE
𝒩 = Helmholtz2D.hypersingular(; wavenumber=k) # Helmholtz2D.hypersingular(;alpha=im*k*η0, beta=+im/(k*η0), wavenumber=k)
M_TEEFIE = assemble(𝒩, X1, X1)


# 1. Excitation: Planewave
H0 = 1.0 
Hz_pw_inc = Helmholtz2D.planewave(; amplitude=H0, wavenumber=k, direction=SVector(1.0, 0.0))
hz_pw_inc = -assemble(DirichletTrace(Hz_pw_inc), X1)
j_TEMFIE_pw = M_TEMFIE \ hz_pw_inc


Hz_pw_sca_num = -potential(HH2DDoubleLayerNear(im * k), pts, j_TEMFIE_pw, X1; type=ComplexF64)
Hz_pw_sca_ana = TE_pec_planewave_H(H0, k, a, pts)

norm(Hz_pw_sca_num - Hz_pw_sca_ana) /norm(Hz_pw_sca_ana)

Et_pw_inc = 1 / (im * ω * ε0) * curl(Hz_pw_inc)
@test Eyinc.direction == Hzinc.direction
@test Eyinc.polarization == SVector(0.0, 1.0)
@test Eyinc.amplitude ≈ Eyinc.gamma / (im * ω * ε0) atol=1e-15
et_pw_inc = assemble(TangentTrace(Et_pw_inc), X1)

j_TEEFIE_pw = (1 / (im * ω * ε0)  * M_TEEFIE) \ et_pw_inc

norm(j_TEEFIE_pw - j_TEMFIE_pw)/norm(j_TEMFIE_pw)

Et_pw_sca_num = 1 / (im * ω * ε0) * potential(HH2DHyperSingularNear(im * k), pts, j_TEEFIE_pw, X1; type=SVector{2, ComplexF64})
Et_pw_sca_ana = TE_pec_planewave_E(H0, k, a, pts)
## We compute the scattered Ez component (scalar)
Et_pw_sca_num[end-3]
Et_pw_sca_ana[end-3]

