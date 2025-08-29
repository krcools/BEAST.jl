using CompScienceMeshes, BEAST


T = tetmeshsphere(1.0,0.3)
X = nedelecc3d(T)
Γ = boundary(T)

X, Y = raviartthomas(Γ), buffachristiansen(Γ)

ω = 1.0
ϵ, ϵ′  = 1.0, 5.0
μ, μ′  = 1.0, 1.0
κ,  η  = ω*√(ϵ*μ), √(μ/ϵ)
κ′, η′ = ω*√(ϵ′*μ′), √(μ′/ϵ′)



N = NCross()

T  = Maxwell3D.singlelayer(wavenumber=κ)
T′ = Maxwell3D.singlelayer(wavenumber=κ′)
K  = Maxwell3D.doublelayer(wavenumber=κ)
K′ = Maxwell3D.doublelayer(wavenumber=κ′)

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
H = -1/(im*κ*η)*curl(E)
e, h = (n × E) × n, (n × H) × n

@hilbertspace j m
@hilbertspace k l

α, α′ = μ/η, μ′/η′
mueller = @discretise(
      (ϵ*η*T-ϵ′*η′*T′)[k,j] -    ((ϵ+ϵ′)/2*N+ϵ*K-ϵ′*K′)[k,m] +
((μ+μ′)/2*N+μ*K-μ′*K′)[l,j] +              (α*T-α′*T′)[l,m] == -ϵ*e[k] - μ*h[l],
j∈X, m∈X, k∈Y, l∈Y)
u,A_mueller,rhs_mueller = solve(mueller)

#postprocessing
using Plots
import Plotly
using LinearAlgebra
fcrj, _ = facecurrents(u[j],X)
fcrm, _ = facecurrents(u[m],X)
Plotly.plot(patch(Γ, norm.(fcrj)),Plotly.Layout(title="j Mueller"))
Plotly.plot(patch(Γ, norm.(fcrm)),Plotly.Layout(title="m Mueller"))


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

import Base.Threads: @spawn
task1 = @spawn nearfield(u[m],u[j],X,X,κ,η,nfpoints,E,H)
task2 = @spawn nearfield(-u[m],-u[j],X,X,κ′,η′,nfpoints)

E_ex, H_ex = fetch(task1)
E_in, H_in = fetch(task2)

E_tot = E_in + E_ex
H_tot = H_in + H_ex

contour(real.(getindex.(E_tot,1)))
contour(real.(getindex.(H_tot,2)))

heatmap(Z, Y, real.(getindex.(E_tot,1)),title="E_tot Mueller")
#heatmap(Z, Y, real.(getindex.(H_tot,2)))
#heatmap(Z, Y, real.(getindex.(E_in,1)))
#heatmap(Z, Y, real.(getindex.(E_ex,1)))