using CompScienceMeshes, BEAST
using LinearAlgebra

Γ = meshcuboid(1.0, 1.0, 1.0, 0.15)
X = raviartthomas(Γ)

κ,  η  = 1.0, 1.0
κ′, η′ = 1.4κ, η/1.4

T  = Maxwell3D.singlelayer(wavenumber=κ)
T′ = Maxwell3D.singlelayer(wavenumber=κ′)
K  = Maxwell3D.doublelayer(wavenumber=κ)
K′ = Maxwell3D.doublelayer(wavenumber=κ′)

d = normalize(x̂+ŷ+ẑ)
p = normalize(d × ẑ)
E = Maxwell3D.planewave(direction=d, polarization=p, wavenumber=κ)
H = -1/(im*κ*η)*curl(E)

e, h = (n × E) × n, (n × H) × n

@hilbertspace j m
@hilbertspace k l

α, α′ = 1/η, 1/η′
pmchwt = @discretise(
    (η*T+η′*T′)[k,j] +      (K+K′)[k,m] -
         (K+K′)[l,j] + (α*T+α′*T′)[l,m] == h[k] + e[l],
    j∈X, m∈X, k∈X, l∈X)

u = solve(pmchwt)

Θ, Φ = range(0.0,stop=2π,length=100), 0.0
ffpoints = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]
farfield = potential(MWFarField3D(κ*im), ffpoints, u, X)

using Plots
using LinearAlgebra
plot(Θ,norm.(farfield))
