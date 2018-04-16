using BEAST, CompScienceMeshes, MieSeries

Θ = linspace(0, π, 100)
Φ = linspace(0, 0,   1)
P = [ [cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)] for θ in Θ for ϕ in Φ]

## Solve the scatering problem using BEM
Γ = readmesh(joinpath(dirname(@__FILE__),"sphere2.in"))
X = raviartthomas(Γ)

ω = κ = 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

@hilbertspace j; @hilbertspace k
efie = @discretise t[k,j]==e[k]  j∈X k∈X
u = gmres(efie)
farbem = 1/(4pi)* norm.(potential(MWFarField3D(im*κ), P, u, X))

## Solve the scattering problem by computing the Mie series
sphere = PECSphere(1.0, 1.0, 1.0)
exc = excite(sphere, ω, :x)
farmie = [norm(collect(farfield(sphere, exc, p)[1:3])) for p in P]

# Compare
@show norm(farmie - farbem) / norm(farmie)
