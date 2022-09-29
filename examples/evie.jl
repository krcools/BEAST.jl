using CompScienceMeshes, BEAST
using LinearAlgebra
using Profile
using StaticArrays

function tau(x::SVector{U,T}) where {U,T}
    4.0-1.0
end

ttrc = X->ttrace(X,y)

T = CompScienceMeshes.tetmeshsphere(1.0,0.25)
X = nedelecc3d(T)
y = boundary(T)
@show numfunctions(X)

ϵ, μ, ω = 1.0, 1.0, 1.0; κ, η = ω * √(ϵ*μ), √(μ/ϵ)
ϵ_r = 4.0
χ = tau
K, I, B = VIE.singlelayer2(wavenumber=κ, tau=χ), Identity(), VIE.boundary2(wavenumber=κ, tau=χ)
E = VIE.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
H = -1/(im*κ*η)*curl(E)

@hilbertspace j 
@hilbertspace k
α = ϵ_r
eq = @varform α*I[j,k]-K[j,k]+B[ttrc(j),k] == E[j]  

evie = @discretise eq j∈X k∈X
u = solve(evie)

#postprocessing
Φ, Θ = [0.0], range(0,stop=2π,length=100)
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ ]
ffe = potential(VIE.farfield(wavenumber=κ, tau=χ), pts, u, X)
ff = ffe

using Plots

#Farfield
plot(xlabel="theta")
plot!(Θ, norm.(ff), label="far field", title="E-VIE")


#Nearfield
Z = range(-1,1,length=100)
Y = range(-1,1,length=100)
nfpoints = [point(0,y,z) for  y in Y, z in Z]

Enear = BEAST.grideval(nfpoints,u,X)
Enear = reshape(Enear,100,100)

contour(real.(getindex.(Enear,1)))
heatmap(Z, Y, real.(getindex.(Enear,1)), clim=(0.0,1.0))

