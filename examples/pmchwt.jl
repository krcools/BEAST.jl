using CompScienceMeshes, BEAST
using LinearAlgebra

Γ = meshcuboid(0.5, 1.0, 1.0, 0.075)
CompScienceMeshes.translate!(Γ, point(-0.25, -0.5, -0.5))
# Γ = meshsphere(1.0, 0.075)
X = raviartthomas(Γ)
@show numfunctions(X)

κ,  η  = 6.0, 1.0
κ′, η′ = 2.4κ, η/2.4

# κ = 1.0
# η = 1.0
# κ′ = 1.5*κ
# η′ = η/1.5

T  = Maxwell3D.singlelayer(wavenumber=κ)
T′ = Maxwell3D.singlelayer(wavenumber=κ′)
K  = Maxwell3D.doublelayer(wavenumber=κ)
K′ = Maxwell3D.doublelayer(wavenumber=κ′)

# d = normalize(x̂+ŷ+ẑ)
d = ẑ
# p = normalize(d × ẑ)
p = x̂
E = Maxwell3D.planewave(direction=d, polarization=p, wavenumber=κ)
H = -1/(im*κ*η)*curl(E)

e, h = (n × E) × n, (n × H) × n

@hilbertspace j m
@hilbertspace k l

α, α′ = 1/η, 1/η′
pmchwt = @discretise(
    (η*T+η′*T′)[k,j] -      (K+K′)[k,m] +
         (K+K′)[l,j] + (α*T+α′*T′)[l,m] == -e[k] - h[l],
    j∈X, m∈X, k∈X, l∈X)

u = solve(pmchwt)

nX = numfunctions(X)
uj = u[1:nX]
um = u[nX+1:end]

Θ, Φ = range(0.0,stop=2π,length=100), 0.0
ffpoints = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]

# Don't forget the far field comprises two contributions
ffm = potential(MWFarField3D(κ*im), ffpoints, um, X)
ffj = potential(MWFarField3D(κ*im), ffpoints, uj, X)
ff = η*im*κ*ffj + im*κ*cross.(ffpoints, ffj)

using Plots
plot(xlabel="theta")
plot!(Θ,norm.(ff),label="far field")

import Plotly
using LinearAlgebra
fcrj, _ = facecurrents(uj,X)
fcrm, _ = facecurrents(um,X)
Plotly.plot(patch(Γ, norm.(fcrj)))
Plotly.plot(patch(Γ, norm.(fcrm)))


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

Z = range(-1,1,length=100)
Y = range(-1,1,length=100)
nfpoints = [point(0,y,z) for z in Z, y in Y]

# nfm_in = potential(BEAST.MWDoubleLayerField3D(wavenumber=κ′), nfpoints, um, X)
# nfj_in = potential(BEAST.MWSingleLayerField3D(wavenumber=κ′), nfpoints, uj, X)
# E_in = nfm_in - η′ * nfj_in

# nfm_ex = potential(BEAST.MWDoubleLayerField3D(wavenumber=κ), nfpoints, um, X)
# nfj_ex = potential(BEAST.MWSingleLayerField3D(wavenumber=κ), nfpoints, uj, X)
# E_ex = -nfm_ex + η * nfj_ex + E.(nfpoints)

E_ex, H_ex = nearfield(um,uj,X,X,κ,η,nfpoints,E,H)
E_in, H_in = nearfield(-um,-uj,X,X,κ′,η′,nfpoints)

E_tot = E_in + E_ex
H_tot = H_in + H_ex

contour(real.(getindex.(E_tot,1)))
heatmap(Z, Y, real.(getindex.(E_tot,1)))
plot(real.(getindex.(E_tot[:,51],1)))



# Compute also the near magnetic field
# Hm_ex = 1/η*potential(BEAST.MWSingleLayerField3D(wavenumber=κ), nfpoints, um, X)
# Hj_ex = potential(BEAST.MWDoubleLayerField3D(wavenumber=κ), nfpoints, uj, X)
# H_ex = Hm_ex + Hj_ex + H.(nfpoints)

# Hm_in = -1/η′*potential(BEAST.MWSingleLayerField3D(wavenumber=κ′), nfpoints, um, X)
# Hj_in = -potential(BEAST.MWDoubleLayerField3D(wavenumber=κ′), nfpoints, uj, X)
# H_in = Hm_in + Hj_in
# Hnf = H_ex + H_in 

contour(real.(getindex.(H_tot,2)))
heatmap(Z, Y, real.(getindex.(H_tot,2)))
plot!(real.(getindex.(H_tot[:,51],2)))