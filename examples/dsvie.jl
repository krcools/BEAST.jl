using CompScienceMeshes, BEAST
using LinearAlgebra
using Profile
using StaticArrays



ntrc = X->ntrace(X,Γ)

T = tetmeshsphere(1.0,0.5)
X = nedelecd3d(T)
Γ = boundary(T)
Y = raviartthomas(Γ)

@show numfunctions(X)
@show numfunctions(Y)


κ,  η  = 1.0, 1.0
κ′, η′ = √2.0κ, η/√2.0
ϵ_r =2.0


χ = x->(1.0-1.0/ϵ_r)

#Volume-Volume
L,I,B = VIE.singlelayer(wavenumber=κ', tau=χ), Identity(), VIE.boundary(wavenumber=κ', tau=χ)
#Volume-Surface 
Lt,Bt,Kt = transpose(VSIE.singlelayer(wavenumber=κ', tau=χ)), transpose(VSIE.boundary(wavenumber=κ', tau=χ)), transpose(VSIE.doublelayer(wavenumber=κ', tau=χ))
Ls,Bs,Ks = VSIE.singlelayer(wavenumber=κ'), VSIE.boundary(wavenumber=κ'), VSIE.doublelayer(wavenumber=κ')
#Surface-Surface
T  = Maxwell3D.singlelayer(wavenumber=κ)  #Outside
T′ = Maxwell3D.singlelayer(wavenumber=κ′) #Inside
K  = Maxwell3D.doublelayer(wavenumber=κ)  #Outside
K′ = Maxwell3D.doublelayer(wavenumber=κ′) #Inside

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
H = -1/(im*κ*η)*curl(E)

e, h = (n × E) × n, (n × H) × n

@hilbertspace D j m
@hilbertspace k l o
β = 1.0/ϵ_r
α, α′ = 1/η, 1/η′

eq = @varform (β*I[k,D]-L[k,D]-B[ntrc(k),D] +  Lt[k,j]+Bt[ntrc(k),j] +          Kt[k,m] +
              Ls[l,D]+Bs[l,ntrc(D)] +       (η*T+η′*T′)[l,j] -      (K+K′)[l,m] +
                         Ks[o,D] +            (K+K′)[o,j] + (α*T+α′*T′)[o,m] == -e[l] - h[o])
  

dvsie = @discretise eq  D∈X k∈X j∈Y m∈Y l∈Y o∈Y

u = solve(dvsie)


#Post processing
Θ, Φ = range(0.0,stop=2π,length=100), 0.0
ffpoints = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]

# Don't forget the far field comprises two contributions
ffm = potential(MWFarField3D(κ*im), ffpoints, u[m], Y)
ffj = potential(MWFarField3D(κ*im), ffpoints, u[j], Y)
ff = -η*im*κ*ffj + im*κ*cross.(ffpoints, ffm)


ffd = potential(VIE.farfield(wavenumber=κ, tau=χ), ffpoints, u[D], X)
ff2 = im*κ*ffd

using Plots
plot(xlabel="theta")
plot!(Θ,norm.(ff),label="far field",title="DVSIE")
plot!(Θ,norm.(ff2),label="far field",title="DVSIE 2")
