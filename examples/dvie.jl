using CompScienceMeshes, BEAST

using LinearAlgebra
#using Profile
#using TimerOutputs
using StaticArrays
using Plots

function tau(x::SVector{U,T}) where {U,T}
    1.0-1.0/5.0
end

ntrc = X->ntrace(X,y)

T = tetmeshsphere(1.0,0.2)
X = nedelecd3d(T)
y = boundary(T)
@show numfunctions(X)

ϵ, μ, ω = 1.0, 1.0, 1.0; κ, η = ω * √(ϵ*μ), √(μ/ϵ)
ϵ_r =5.0
χ = tau
K, I, B = VIE.singlelayer(wavenumber=κ, tau=χ), Identity(), VIE.boundary(wavenumber=κ, tau=χ)
E = VIE.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
H = -1/(im*κ*η)*curl(E)

@hilbertspace j 
@hilbertspace k
α = 1.0/ϵ_r
eq = @varform α*I[j,k]-K[j,k]-B[ntrc(j),k] == E[j] 

dbvie = @discretise eq j∈X k∈X
u_m = solve(dbvie)
#=
#postprocessing
Φ, Θ = [0.0], range(0,stop=2π,length=100)
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]
ffd = potential(VIE.farfield(wavenumber=κ, tau=χ), pts, u, X)
ff = im*κ*ffd

using Plots

#Farfield
plot(xlabel="theta")
plot!(Θ, norm.(ff), label="far field", title="D-VIE")
=#

#NearField
Z = range(-1,1,length=100)
Y = range(-1,1,length=100)
nfpoints = [point(0,y,z) for z in Z, y in Y]

Enear = BEAST.grideval(nfpoints,α.*u_m,X)
Enear = reshape(Enear,100,100)

contour(real.(getindex.(Enear,1)))
heatmap(Z, Y,  real.(getindex.(Enear,1)))


plot!(Y[2:99],real.(getindex.(Enear[2:99,50],1)),label="D-VIE (simplex)", linestyle=:dash, linecolor=:darkorange4)
#plot!(Y[2:99],real.(getindex.(Enear_simplex[2:99,50],1)),label="D-VIE (simplex)", linestyle=:solid, linecolor=:darkorange3)

#=
import Cairo
using DataFrames, Gadfly, RDatasets
D = dataset("datasets","HairEyeColor")
palette = ["skyblue","skyblue3"]


D = DataFrame(
   Type=["D-VIE \n (1,190 Tets.)","D-VIE \n (1,190 Tets.)","D-VIE \n (1,190 Tets.)","D-VIE \n (1,190 Tets.)","E-VIE \n (1,190 Tets.)","E-VIE \n (1,190 Tets.)","E-VIE \n (1,190 Tets.)","E-VIE \n (1,190 Tets.)","D-VIE \n (2,706 Tets.)","D-VIE \n (2,706 Tets.)","D-VIE \n (2,706 Tets.)","D-VIE \n (2,706 Tets.)","E-VIE \n (2,706 Tets.)","E-VIE \n (2,706 Tets.)","E-VIE \n (2,706 Tets.)","E-VIE \n (2,706 Tets.)"], 
   Quad=["Simplex","Simplex","Classic","Classic","Simplex","Simplex","Classic","Classic","Simplex","Simplex","Classic","Classic","Simplex","Simplex","Classic","Classic"],
   Sing=["Regular","Singular","Regular","Singular","Regular","Singular","Regular","Singular","Regular","Singular","Regular","Singular","Regular","Singular","Regular","Singular"],
   Time=[24,4.77,23.1,27.9,41.9,8.49,41.5,45.7,125,13.1,125,68.6,222,20.5,222,115],
)

p = Gadfly.plot(D, x=:Quad, y=:Time, color=:Sing, xgroup=:Type, 
    Geom.subplot_grid(Geom.bar(position=:stack)),
    Scale.color_discrete_manual(palette...), 
    Guide.xlabel(""),Guide.ylabel("Assembly time [s]"), Guide.colorkey("Interaction"),Theme(major_label_font="Times",minor_label_font="Times",key_title_font="Times", key_label_font="Times",background_color=color("white")))

# draw(PNG("haireyecolor.png", 6.6inch, 4inch), p)
draw(PNG(6.6inch, 4inch), p)
=#