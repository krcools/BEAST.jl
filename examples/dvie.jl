using CompScienceMeshes, BEAST
using LinearAlgebra
using Profile
using StaticArrays

function tau(x::SVector{U,T}) where {U,T}
    1.0-1.0/4.0
end

ntrc = X->ntrace(X,y)

T = CompScienceMeshes.tetmeshsphere(1.0,0.25)
X = nedelecd3d(T)
y = boundary(T)
@show numfunctions(X)

ϵ, μ, ω = 1.0, 1.0, 1.0; κ, η = ω * √(ϵ*μ), √(μ/ϵ)
ϵ_r =4.0
χ = tau
K, I, B = VIE.singlelayer(wavenumber=κ, tau=χ), Identity(), VIE.boundary(wavenumber=κ, tau=χ)
E = VIE.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
H = -1/(im*κ*η)*curl(E)

@hilbertspace j 
@hilbertspace k
α = 1.0/ϵ_r
eq = @varform α*I[j,k]-K[j,k]-B[ntrc(j),k] == E[j] 

dbvie = @discretise eq j∈X k∈X
u = solve(dbvie)

#postprocessing
Φ, Θ = [0.0], range(0,stop=2π,length=100)
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]
ffd = potential(VIE.farfield(wavenumber=κ, tau=χ), pts, u, X)
ff = im*κ*ffd

using Plots

#Farfield
plot(xlabel="theta")
plot!(Θ, norm.(ff), label="far field", title="D-VIE")

#NearField
Z = range(-1,1,length=100)
Y = range(-1,1,length=100)
nfpoints = [point(0,y,z) for  y in Y, z in Z]

Enear = BEAST.grideval(nfpoints,α.*u,X)
Enear = reshape(Enear,100,100)

contour(real.(getindex.(Enear,1)))
heatmap(Z, Y,  real.(getindex.(Enear,1)), clim=(0.0,1.0))




