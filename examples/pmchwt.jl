using CompScienceMeshes, BEAST
using LinearAlgebra

# Γ = meshcuboid(1.0, 1.0, 1.0, 0.15)
Γ = meshsphere(1.0, 0.075)
X = raviartthomas(Γ)
@show numfunctions(X)

κ,  η  = 6.0, 1.0
κ′, η′ = 2.4κ, η/2.4

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
         (K+K′)[l,j] + (α*T+α′*T′)[l,m] == e[k] - h[l],
    j∈X, m∈X, k∈X, l∈X)

u = solve(pmchwt)

nX = numfunctions(X)
uj = u[1:nX]
um = u[nX+1:end]

Θ, Φ = range(0.0,stop=2π,length=100), 0.0
ffpoints = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]

# Don't forgt the far field comprises two contributions
ffm = potential(MWFarField3D(κ*im), ffpoints, um, X)
ffj = potential(MWFarField3D(κ*im), ffpoints, uj, X)
ff = η*im*κ*ffj + im*κ*cross.(ffpoints, ffm)

using Plots
plot(xlabel="theta")
plot!(Θ,norm.(ff),label="far field")

import Plotly
using LinearAlgebra
fcrj, _ = facecurrents(uj,X)
fcrm, _ = facecurrents(um,X)
Plotly.plot(patch(Γ, norm.(fcrm)))

Z = range(-2,2,length=100)
Y = range(-2,2,length=100)
nfpoints = [point(0,y,z) for z in Z, y in Y]
nfm_in = potential(BEAST.MWDoubleLayerField3D(wavenumber=κ′), nfpoints, um, X)
nfj_in = potential(BEAST.MWSingleLayerField3D(wavenumber=κ′), nfpoints, uj, X)
nf_in = -nfm_in + η′ * nfj_in

nfm_ex = potential(BEAST.MWDoubleLayerField3D(wavenumber=κ), nfpoints, um, X)
nfj_ex = potential(BEAST.MWSingleLayerField3D(wavenumber=κ), nfpoints, uj, X)
nf_ex = nfm_ex - η * nfj_ex

nf_in = reshape(nf_in, size(nfpoints))
nf_ex = reshape(nf_ex, size(nfpoints))
nf = nf_in + nf_ex + E.(nfpoints)

contour(real.(getindex.(nf,1)))

plot()
plot(real.(getindex.(nf[:,51],1)))
