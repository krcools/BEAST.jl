using Test

using CompScienceMeshes
using BEAST
using StaticArrays
using LinearAlgebra

c = 3e8
μ = 4*π*1e-7
ε = 1/(μ*c^2)
f = 5e7
λ = c/f
k = 2*π/λ
ω = k*c
η = sqrt(μ/ε)

a = 1
Γ_orig = CompScienceMeshes.meshcuboid(a,a,a,0.1)
Γ = translate(Γ_orig,SVector(-a/2,-a/2,-a/2))

Φ, Θ = [0.0], range(0,stop=π,length=100)
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]

# This is an electric dipole
# The pre-factor (1/ε) is used to resemble 
# (9.18) in Jackson's Classical Electrodynamics
E = (1/ε) * dipolemw3d(location=SVector(0.4,0.2,0), 
                       orientation=1e-9.*SVector(0.5,0.5,0), 
                       wavenumber=k)

n = BEAST.NormalVector()

𝒆 = (n × E) × n
H = (-1/(im*μ*ω))*curl(E)
𝒉 = (n × H) × n

𝓣 = Maxwell3D.singlelayer(wavenumber=k)

X = raviartthomas(Γ)

T = Matrix(assemble(𝓣,X,X))
e = Vector(assemble(𝒆,X))
j_EFIE = T\e

nf_E_EFIE = potential(MWSingleLayerField3D(wavenumber=k), pts, j_EFIE, X)
nf_H_EFIE = potential(BEAST.MWDoubleLayerField3D(wavenumber=k), pts, j_EFIE, X) ./ η
ff_E_EFIE = potential(MWFarField3D(wavenumber=k), pts, j_EFIE, X)

accuracy_default_quadrature_nf_e = norm(nf_E_EFIE - E.(pts))/norm(E.(pts))
accuracy_default_quadrature_nf_h = norm(nf_H_EFIE - H.(pts))/norm(H.(pts))
accuracy_default_quadrature_ff_e = 
    norm(ff_E_EFIE - E.(pts, isfarfield=true))/norm(E.(pts, isfarfield=true))

function fastquaddata(op::BEAST.MaxwellOperator3D,
    test_local_space::BEAST.RefSpace, trial_local_space::BEAST.RefSpace,
    test_charts, trial_charts)

    a, b = 0.0, 1.0
    # CommonVertex, CommonEdge, CommonFace rules
    println("Fast quadrule is called")
    tqd = quadpoints(test_local_space, test_charts, (1,2))
    bqd = quadpoints(trial_local_space, trial_charts, (1,2))
    leg = (BEAST._legendre(3,a,b), BEAST._legendre(4,a,b), BEAST._legendre(5,a,b),)

    return (tpoints=tqd, bpoints=bqd, gausslegendre=leg)
end

function fastquaddata(fn::BEAST.Functional, refs, cells) 
    println("Fast RHS quadrature")
    return quadpoints(refs, cells, [1])
end


fastT = Matrix(assemble(𝓣,X,X, quaddata=fastquaddata))
faste = Vector(assemble(𝒆,X, quaddata=fastquaddata))
fastj_EFIE = fastT\faste

fastnf_E_EFIE = potential(MWSingleLayerField3D(wavenumber=k), pts, fastj_EFIE, X)
fastnf_H_EFIE = potential(BEAST.MWDoubleLayerField3D(wavenumber=k), pts, fastj_EFIE, X) ./ η
fastff_E_EFIE = potential(MWFarField3D(wavenumber=k), pts, fastj_EFIE, X)

accuracy_fast_quadrature_nf_e = norm(fastnf_E_EFIE - E.(pts))/norm(E.(pts))
accuracy_fast_quadrature_nf_h = norm(fastnf_H_EFIE - H.(pts))/norm(H.(pts))
accuracy_fast_quadrature_ff_e = 
    norm(fastff_E_EFIE - E.(pts, isfarfield=true))/norm(E.(pts, isfarfield=true))

@test accuracy_fast_quadrature_nf_e > accuracy_default_quadrature_nf_e
@test accuracy_fast_quadrature_nf_h > accuracy_default_quadrature_nf_h
@test accuracy_fast_quadrature_ff_e > accuracy_default_quadrature_ff_e

