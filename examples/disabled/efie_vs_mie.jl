using BEAST, CompScienceMeshes, MieSeries

Θ = range(0, stop=π, length=100)
Φ = 0.0
P = [ [cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)] for θ in Θ for ϕ in Φ]

## Solve the scatering problem using BEM
Γ = readmesh(joinpath(dirname(@__FILE__),"sphere2.in"))
X = raviartthomas(Γ)

ϵ = μ = 1.0
ω = 1.0; κ = ω * √(ϵ*μ)
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

@hilbertspace j; @hilbertspace k
efie = @discretise t[k,j]==e[k]  j∈X k∈X
u = gmres(efie)

A = Maxwell3D.farfield(wavenumber=κ)
farbem = 1/(4pi)* norm.(potential(A, P, u, X))

## Solve the scattering problem by computing the Mie series
radius = 1.0; ϵ = μ = 1.0
sphere = PECSphere(radius, ϵ, μ)
exc = excite(sphere, ω, :x)
farmie = [norm(collect(farfield(sphere, exc, p)[1:3])) for p in P]

# Compare
@show norm(farmie - farbem) / norm(farmie)
