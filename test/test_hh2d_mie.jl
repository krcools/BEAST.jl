using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test
using SpecialFunctions

const ε0 = 8.854187821e-12
const μ0 = 4π*1e-7
const c0 = 1/sqrt(ε0*μ0)
const η = sqrt(μ0/ε0)

f = 1e9 # 1 GHz
ω = 2π*f
λ = c0/f
k = 2π/λ
h = λ/10
ε
μ
a = 1.0 # radius of the scatterer

circle = CompScienceMeshes.meshcircle(a, h)
X0 = lagrangecxd0(circle)

##############
# a) TM-EFIE

S = Helmholtz2D.singlelayer(;alpha=im*k*η, wavenumber=k)
𝗦 = assemble(S,X0,X0)


##
# 1. RHS is a planewave

#=
 We consider the case of scattering by a circular PEC cylinder [Sec 6.4.2, Jin,
 Theory and Computation of Electromagnetic Fields].

 The incoming planewave is along the x-axis.
=#
E0 = 1.0 # amplitude
Einc = Helmholtz2D.planewave(; amplitude=E0, wavenumber=k, direction=SVector(1.0, 0.0))

𝗲 = assemble(DirichletTrace(Einc), X0)
𝗷 = 𝗦 \ 𝗲

# Computing the fields on a circle of radius rc
rc = 2.0*a
pts = meshcircle(rc, 0.6 * rc).vertices

# We compute the scattered Ez component (scalar)
Esca_num = -im*k*η * potential(HH2DSingleLayerNear(im * k), pts, 𝗷, X0; type=ComplexF64)

#=
 For the analytical evaluation, we consider [Eq. 6.4.12, Jin] as the expression for the
 scattered Electric field.
=#
function TME_pec_planewave(E0,k,a,ρ,φ)
    n = 0
    val(n) = (1.0im)^(-n) * (besselj(n,k*a)/hankelh2(n,k*a))*hankelh2(n,k*ρ)*exp(im*n*φ)
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
    return -E0 * ret
end

function TME_pec_planewave(E0,k,a,pts)
    cart2polar(x,y) = SVector(sqrt(x^2+y^2), atan(y,x))
    return [TME_pec_planewave(E0,k,a,cart2polar(p[1],p[2])...) for p in pts]
end

Esca_ana = TME_pec_planewave(E0,k,a,pts)

abs.(Esca_num - Esca_ana) ./ abs.(Esca_ana)
norm(Esca_num - Esca_ana) ./ norm(Esca_ana)

##
# 2. RHS is an infinite line current

#=
 Here, the line current is infinitely long along the z-axis, and is pointed towards the 
 z-axis as well. Thus, the potentials Ax, Ay will not be excited, leaving us only Az to 
 compute, which will then give us Ez. The monopole excitation in BEAST is used for this 
 excitation.
=#

# Ref: Sec 6.5.3, Jin, Theory and Computation of Electromagnetic fields
# Analytical expression used for scattered field: Eq 6.5.38

# Analytical result
function TME_pec_line_curr(I,k,a,ρ,φ,ρp,φp)
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
    return η*k*I/4 * ret
end

function TME_pec_line_curr(I,k,a,pts,spt)
    cart2polar(x,y) = SVector(sqrt(x^2+y^2), atan(y,x))
    return [TME_pec_line_curr(I,k,a,cart2polar(p[1],p[2])..., spt[1],spt[2]) for p in pts]
end
# Rhs: Line current: Monopole
# Numerical evaluation

I = 1.0 # 1 A
ρp = 2.0 # Source radial pos
φp = 0.0 # Source angular pos
xp = ρp * cos(φp)
yp = ρp * sin(φp)

# Choosing the amplitude of Einc such that [Eq 6.5.11, Jin] is satisfied
Einc = Helmholtz2D.monopole(;position = SVector(xp,yp), wavenumber=k, amplitude = -η*k*I/4)
𝗲 = assemble(DirichletTrace(Einc),X0)
𝗷 = 𝗦 \ 𝗲

pts = meshcircle(3.0*a, 3.0*0.6*a).vertices
Esca_num = -im*η*k * potential(HH2DSingleLayerNear(im*k),pts,𝗷,X0;type=ComplexF64)

Esca_ana = TME_pec_line_curr(I,k,a,pts,SVector(ρp,φp))

abs.(Esca_num - Esca_ana) ./ abs.(Esca_ana)
norm(Esca_num - Esca_ana) / norm(Esca_ana)



## TE-MFIE and TE-EFIE

H0 = 1.0
 
Hzinc = Helmholtz2D.planewave(; amplitude=H0, wavenumber=k, direction=SVector(1.0, 0.0))

Eyinc = 1 / (im * ω * ε0) * curl(Hzinc)

@test Eyinc.direction == Hzinc.direction
@test Eyinc.polarization == SVector(0.0, 1.0)
@test Eyinc.amplitude ≈ Eyinc.gamma / (im * ω * ε0) atol=1e-15

X1 = lagrangec0d1(circle)


# TE-MFIE
#0.5I - D

D = Helmholtz2D.doublelayer(; wavenumber=k) #Helmholtz2D.hypersingular(;alpha=im*k*η, beta=-im/(k*η), wavenumber=k)
𝗵 = -assemble(DirichletTrace(Hzinc), X1)
𝗜 = assemble(BEAST.Identity(), X1, X1)
𝗗 = assemble(D, X1, X1)

jTEMFIE = (0.5*𝗜 - 𝗗)\𝗵

Hzsca_num = potential(HH2DDoubleLayerNear(im * k), pts, jTEMFIE, X1; type=ComplexF64)

# TE-EFIE
N = Helmholtz2D.hypersingular(; wavenumber=k) #Helmholtz2D.hypersingular(;alpha=im*k*η, beta=-im/(k*η), wavenumber=k)
𝗡 = assemble(N, X1, X1)
𝗲t = (im * ω * ε0) * assemble(TangentTrace(Eyinc), X1)
jTEEFIE = 𝗡 \ 𝗲t

# We compute the scattered Ez component (scalar)
Esca_num = -1 / (im * ω * ε0) * potential(HH2DHyperSingularNear(im * k), pts, 𝗷t, X1; type=SVector{2,ComplexF64})

dbesselj(n,x) = besselj(n-1,x) - n/x * besselj(n,x)
dhankelh2(n,x) = hankelh2(n-1,x) - n/x * hankelh2(n,x)

function TEE_pec_planewave(H0, k, a, ρ, φ)
    n = 0

    # Jin (6.4.19)
    b_p(n) = -(1.0im)^(-n)*dbesselj(n,k*a)/dhankelh2(n,k*a)

    # Hz is not needed, but maybe we use it later for MFIE
    valHz(n) = b_p(n) * hankelh2(n, k*ρ)*exp(im*n*φ)
    valEφ(n) = -1/(im * ω * ε0) * k * b_p(n) * dhankelh2(n, k*ρ)*exp(im*n*φ)
    valEρ(n) = 1/(im * ω * ε0)  * (im * n / ρ) * b_p(n) * hankelh2(n, k*ρ)*exp(im*n*φ)

    retHz = valHz(0)
    retEφ = valEφ(0)
    retEρ = valEρ(0)
    while true
        n += 1
        bufHz = valHz(n) + valHz(-n)
        bufEφ = valEφ(n) + valEφ(-n)
        bufEρ = valEρ(n) + valEρ(-n)
        if abs(bufHz) >= abs(retHz) * eps(eltype(real(retHz))) * 1e3 && 
            abs(bufEφ) >= abs(retEφ) * eps(eltype(real(retEφ))) * 1e3  && 
            abs(bufEρ) >= abs(retEρ) * eps(eltype(real(retEρ))) * 1e3

            retHz += bufHz
            retEφ += bufEφ
            retEρ += bufEρ
        else
            break
        end
    end

    return retHz, H0*(retEρ*SVector(cos(φ), sin(φ)) + retEφ*SVector(-sin(φ), cos(φ)))
end

function TEE_pec_planewave_H(H0, k, a, pts)
    cart2polar(x,y) = SVector(sqrt(x^2 + y^2), atan(y, x))
    return [TEE_pec_planewave(H0, k, a, cart2polar(p[1],p[2])...) for p in pts]
end

function TEE_pec_planewave_E(H0, k, a, pts)
    cart2polar(x,y) = SVector(sqrt(x^2 + y^2), atan(y, x))
    return [TEE_pec_planewave(H0, k, a, cart2polar(p[1],p[2])...) for p in pts]
end

Esca_ana = TEE_pec_planewave(H0, k, a, pts)

Esca_ana