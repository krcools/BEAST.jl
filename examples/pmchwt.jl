using CompScienceMeshes, BEAST
using LinearAlgebra

T = tetmeshsphere(1.0,0.25)
X = nedelecc3d(T)
Γ = boundary(T)

X = raviartthomas(Γ)
@show numfunctions(X)

κ,  η  = 1.0, 1.0
κ′, η′ = 2.0κ, η/2.0

T  = Maxwell3D.singlelayer(wavenumber=κ)
T′ = Maxwell3D.singlelayer(wavenumber=κ′)
K  = Maxwell3D.doublelayer(wavenumber=κ)
K′ = Maxwell3D.doublelayer(wavenumber=κ′)

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
H = -1/(im*κ*η)*curl(E)

e, h = (n × E) × n, (n × H) × n

@hilbertspace j m
@hilbertspace k l

α, α′ = 1/η, 1/η′
pmchwt = @discretise(
    (η*T+η′*T′)[k,j] -      (K+K′)[k,m] +
         (K+K′)[l,j] + (α*T+α′*T′)[l,m] == -e[k] - h[l],
    j∈X, m∈X, k∈X, l∈X)

# Z = BEAST.sysmatrix(pmchwt)
# b = BEAST.rhs(pmchwt)
# M = BEAST.convert_to_dense(Z)
# u = M \ b
# error()
# W = BEAST.assemble_hide(pmchwt.equation.lhs, pmchwt.test_space_dict, pmchwt.trial_space_dict)

# using BEAST.BlockArrays
# using BEAST.IterativeSolvers

# u = PseudoBlockVector{ComplexF64}(undef, numfunctions.([X,X]))
# fill!(u,0)

# error()
# u, ch = IterativeSolvers.gmres!(u, Z, b, log=true,  maxiter=1000,
#     restart=1000, reltol=1e-5, verbose=true)

u = solve(pmchwt)

Θ, Φ = range(0.0,stop=2π,length=100), 0.0
ffpoints = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]

# Don't forget the far field comprises two contributions
ffm = potential(MWFarField3D(κ*im), ffpoints, u[m], X)
ffj = potential(MWFarField3D(κ*im), ffpoints, u[j], X)
ff = η*im*κ*ffm + im*κ*cross.(ffpoints, ffj)

using Plots
plot(xlabel="theta")
plot!(Θ,norm.(ff),label="far field",title="PMCHWT")

#import Plotly
#using LinearAlgebra
#fcrj, _ = facecurrents(u[j],X)
#fcrm, _ = facecurrents(u[m],X)
#Plotly.plot(patch(Γ, norm.(fcrj)))
#Plotly.plot(patch(Γ, norm.(fcrm)))


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

E_ex, H_ex = nearfield(u[m],u[j],X,X,κ,η,nfpoints,E,H)
E_in, H_in = nearfield(-u[m],-u[j],X,X,κ′,η′,nfpoints)

E_tot = E_in + E_ex
H_tot = H_in + H_ex

contour(real.(getindex.(E_tot,1)))
contour(real.(getindex.(H_tot,2)))

heatmap(Z, Y, real.(getindex.(E_tot,1)))
#heatmap(Z, Y, real.(getindex.(H_tot,2)))

#plot(real.(getindex.(E_tot[:,51],1)))
#plot!(real.(getindex.(H_tot[:,51],2)))


