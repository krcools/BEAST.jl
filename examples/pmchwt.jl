using CompScienceMeshes, BEAST
using LinearAlgebra


function nearfield(um,uj,Xm,Xj,κ,η,points,
    Einc=(x->point(0,0,0)),
    Hinc=(x->point(0,0,0)))

    K = BEAST.MWDoubleLayerField3D(wavenumber=κ)
    T = BEAST.MWSingleLayerField3D(wavenumber=κ)

    Em = potential(K, points, um, Xm)
    Ej = potential(T, points, uj, Xj)
    E = -Em + η * Ej + Einc.(points)

    Hm = potential(T, points, um, Xm)
    Hj = potential(K, points, uj, Xj)
    H = 1/η*Hm + Hj + Hinc.(points)

    return E, H
end

ϵ0 = 8.854e-12
μ0 = 4π*1e-7
c = 1/√(ϵ0*μ0)

# λ = 2.9979563769321627
# ω = 2π*c/λ
ω = 10.0^8*2*pi

# Ω = CompScienceMeshes.tetmeshsphere(λ,0.1*λ)
Ω = CompScienceMeshes.tetmeshsphere(1.0,0.2)
Γ = boundary(Ω)
X = raviartthomas(Γ)
@show numfunctions(X)

ϵr = 2.0
μr = 3.0

κ, η = ω/c, √(μ0/ϵ0)
κ′, η′ = κ*√(ϵr*μr), η*√(μr/ϵr)

# κ,  η  = π, 1.0
# κ′, η′ = 2.0κ, η/2.0

T  = Maxwell3D.singlelayer(wavenumber=κ)
T′ = Maxwell3D.singlelayer(wavenumber=κ′)
K  = Maxwell3D.doublelayer(wavenumber=κ)
K′ = Maxwell3D.doublelayer(wavenumber=κ′)

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
H = -1/(im*κ*η)*curl(E)

e = (n × E) × n
h = (n × H) × n

@hilbertspace j m
@hilbertspace k l

α, α′ = 1/η, 1/η′
pmchwt = @discretise(
    (η*T+η′*T′)[k,j] -      (K+K′)[k,m] +
         (K+K′)[l,j] + (α*T+α′*T′)[l,m] == -e[k] - h[l],
    j∈X, m∈X, k∈X, l∈X)

u = BEAST.solve(pmchwt)

Θ, Φ = range(0.0,stop=2π,length=100), 0.0
ffpoints = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]

# Don't forget the far field comprises two contributions
ffm = potential(MWFarField3D(κ*im, η), ffpoints, u[m], X)
ffj = potential(MWFarField3D(κ*im, η), ffpoints, u[j], X)                                                                                                                                                                                                                                                          
ff = -η*im*κ*ffj + im*κ*cross.(ffpoints, ffm)

# Compare the far field and the field far
using Plots
ffradius = 100.0
E_far, H_far = nearfield(u[m],u[j],X,X,κ,η, ffradius .* ffpoints)
nxE_far = cross.(ffpoints, E_far) * (4π*ffradius) / exp(-im*κ*ffradius)
Et_far = -cross.(ffpoints, nxE_far)

Plots.plot()
Plots.plot!(Θ, norm.(ff)/η ,label="far field")
Plots.scatter!(Θ, norm.(Et_far), label="field far")

using Plots
Plots.plot(xlabel="theta")
Plots.plot!(Θ,norm.(ff),label="far field",title="PMCHWT")

import Plotly
using LinearAlgebra
fcrj, _ = facecurrents(u[j],X)
fcrm, _ = facecurrents(u[m],X)
Plotly.plot(patch(Γ, norm.(fcrj)))
Plotly.plot(patch(Γ, norm.(fcrm)))





Zz = range(-1,1,length=130)
Yy = range(-1,1,length=100)
nfpoints = [point(z,0,y) for z in Zz, y in Yy]
nfpoints = pts
import Base.Threads: @spawn
task1 = @spawn nearfield(u[m],u[j],X,X,κ,η,nfpoints,E,H)
task2 = @spawn nearfield(-u[m],-u[j],X,X,κ′,η′,nfpoints)

E_ex, H_ex = fetch(task1)
E_in, H_in = fetch(task2)

E_tot = E_in + E_ex
H_tot = H_in + H_ex

Plots.contour(real.(getindex.(E_tot,1)))
Plots.contour(real.(getindex.(H_tot,2)))

Plots.heatmap(Z, Y, clamp.(real.(getindex.(E_tot,1)),-1.5,1.5))
Plots.heatmap(Z, Y, clamp.(imag.(getindex.(E_tot,1)),-1.5,1.5))
Plots.heatmap(Yy, Zz, norm.(getindex.(H_tot,2)))
Plots.heatmap(Y, Z, norm.(getindex.(H_tot,2)))

Plots.plot(real.(getindex.(E_tot[:,51],1)))
Plots.plot(real.(getindex.(H_tot[:,51],2)))

Plots.heatmap(Zz,Yy,norm.(getindex.(-Bbb[1].+Bbb[3].-Bbb[2]/2,2))/μ0,clim=(0,0.02))
Plots.heatmap(Zz,Yy,real(getindex.(H.(pts),2)))
