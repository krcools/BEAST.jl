using CompScienceMeshes, BEAST
using LinearAlgebra
using Profile
#using LiftedMaps
#using BlockArrays
using StaticArrays



ntrc = X->ntrace(X,Γ)

T = tetmeshsphere(1.0,0.3)
X = nedelecd3d(T)
Γ = boundary(T)
Y = raviartthomas(Γ)

@show numfunctions(X)
@show numfunctions(Y)

Z = BEAST.buffachristiansen2(Γ)


κ,  η  = 1.0, 1.0
κ′, η′ = √2κ, η/√2
ϵ_r =3
ϵ_b =2.0
κ_in, η_in = √6κ, η/√6


#χ = x->(1.0-1.0/ϵ_r)

function tau(x::SVector{U,T}) where {U,T}
    1.0-1.0/3.0
end

χ = tau



N = NCross()
#Volume-Volume
L,I,B = VIE.singlelayer(wavenumber=κ′, tau=χ), Identity(), VIE.boundary(wavenumber=κ′, tau=χ)
#Volume-Surface 
Lt,Bt,Kt = transpose(VSIE.singlelayer(wavenumber=κ′)), transpose(VSIE.boundary(wavenumber=κ′)), transpose(VSIE.doublelayer(wavenumber=κ′))
#Kt = VSIE.doublelayerT(wavenumber=κ′)

Ls,Bs,Ks = VSIE.singlelayer(wavenumber=κ′, tau=χ), VSIE.boundary(wavenumber=κ′, tau=χ), VSIE.doublelayer(wavenumber=κ′, tau=χ)
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
β = 1/(ϵ_r*ϵ_b)
ν = 1/ϵ_b
α, α′ = 1/η, 1/η′
γ′ = im*η′/κ′
ζ′ = im*κ′/(ϵ_b*η′)
δ′ = im*κ′/ϵ_b

#=
eq = @varform (β*I[k,D] + η′*Lt[k,j]-γ′*Bt[ntrc(k),j] -          Kt[k,m] +
                                     (η*T+η′*T′)[l,j] -      (K+K′)[l,m] +
                                          (K+K′)[o,j] + (α*T+α′*T′)[o,m] == -e[l] - h[o])

=#

eq = @varform (β*I[k,D]-ν*L[k,D]-ν*B[ntrc(k),D] + η′*Lt[k,j]-γ′*Bt[ntrc(k),j] -          Kt[k,m] +
                    -δ′*Ls[l,D]-ν*Bs[l,ntrc(D)]            + (η*T+η′*T′)[l,j] -      (K+K′)[l,m] +
                                    -ζ′*Ks[o,D] +                 (K+K′)[o,j] + (α*T+α′*T′)[o,m] == -e[l] - h[o])


dvsie = @discretise eq  D∈X k∈X j∈Y m∈Y l∈Y o∈Y

u_n  = solve(dvsie)


#=
#preconditioner
Mxx = assemble(I,X,X)
iMxx = inv(Matrix(Mxx))

Tzz = assemble(T,Z,Z); println("dual discretisation assembled.")
Nyz = Matrix(assemble(N,Y,Z)); println("duality form assembled.")

iNyz = inv(Nyz); println("duality form inverted.")
NTN = iNyz' * Tzz * iNyz 

M = zeros(Int, 3)
N = zeros(Int, 3)

M[1] = numfunctions(X)
N[1] = numfunctions(X)
M[2] = M[3] = numfunctions(Y)
N[2] = N[3] = numfunctions(Y)

U = BlockArrays.blockedrange(M)
V = BlockArrays.blockedrange(N)

precond = BEAST.ZeroMap{Float32}(U, V)

z1 = LiftedMap(iMxx,Block(1),Block(1),U,V)
z2 = LiftedMap(NTN,Block(2),Block(2),U,V)
z3 = LiftedMap(NTN,Block(3),Block(3),U,V)
precond = precond +ϵ_r*z1 + z2 + z3

A_dvsie_precond = precond*A_dvsie

#GMREs
import IterativeSolvers
cT = promote_type(eltype(A_dvsie), eltype(rhs))
x = BlockedVector{cT}(undef, M)
fill!(x, 0)
x, ch = IterativeSolvers.gmres!(x, A_dvsie, Rhs, log=true,  reltol=1e-4)
fill!(x, 0)
x, ch = IterativeSolvers.gmres!(x, precond*A_dvsie, precond*Rhs, log=true,  reltol=1e-8)
=#

#Post processing
Θ, Φ = range(0.0,stop=2π,length=100), 0.0
ffpoints = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]

# Don't forget the far field comprises two contributions
ffm = potential(MWFarField3D(gamma=κ*im), ffpoints, u_n[m], Y)
ffj = potential(MWFarField3D(gamma=κ*im), ffpoints, u_n[j], Y)
ff = -η*im*κ*ffj + im*κ*cross.(ffpoints, ffm)


ffd = potential(VIE.farfield(wavenumber=κ, tau=χ), ffpoints, u_n[D], X)
ff2 = im*κ*ffd

using Plots
plot(xlabel="theta")
plot!(Θ,norm.(ff),label="far field",title="DVSIE")
plot!(Θ,√6.0*norm.(ff),label="far field",title="DVSIE 2")

import Plotly
using LinearAlgebra
fcrj, _ = facecurrents(u_n[j],Y)
fcrm, _ = facecurrents(u_n[m],Y)
vsie_j = Plotly.plot([patch(Γ, norm.(fcrj))], Plotly.Layout(title="j D-VSIE"))
vsie_m = Plotly.plot(patch(Γ, norm.(fcrm)), Plotly.Layout(title="m D-VSIE"))

#NearField
function nearfield(um,uj,Xm,Xj,κ,η,points,
    Einc=(x->point(0,0,0)),
    Hinc=(x->point(0,0,0)))

    K = BEAST.MWDoubleLayerField3D(wavenumber=κ)
    T = BEAST.MWSingleLayerField3D(wavenumber=κ)

    Em = potential(K, points, um, Xm)
    Ej = potential(T, points, uj, Xj)
    E = -Em + η*Ej + Einc.(points)

    Hm = potential(T, points, um, Xm)
    Hj = potential(K, points, uj, Xj)
    H = 1/η*Hm + Hj + Hinc.(points)

    return E, H
end

Zz = range(-1,1,length=100)
Yy = range(-1,1,length=100)
nfpoints = [point(0.0,y,z) for  z in Zz, y in Yy]


import Base.Threads: @spawn
task1 = @spawn nearfield(u_n[m],u_n[j],Y,Y,κ,η,nfpoints,E,H)
task2 = @spawn nearfield(-u_n[m],-u_n[j],Y,Y,κ_in,η_in,nfpoints)

E_ex, H_ex = fetch(task1)
E_in, H_in = fetch(task2)



Enear = BEAST.grideval(nfpoints,β.* u_n[D],X)
Enear = reshape(Enear,100,100)

contour(real.(getindex.(Enear,1)))
heatmap(Zz, Yy,  real.(getindex.(Enear,1)))

heatmap(Zz, Yy, real.(getindex.(E_in,1)))
heatmap(Zz, Yy, real.(getindex.(E_ex,1)))

contour(real.(getindex.(E_ex,1)))


Dd = range(1,100,step=1)
plot!(Yy,real.(getindex.(E_in[Dd,50],1)))
plot(collect(Yy)[2:99],real.(getindex.(Enear[2:99,50],1)))