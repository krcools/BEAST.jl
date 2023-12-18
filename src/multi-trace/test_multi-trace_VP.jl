using BEAST
using CompScienceMeshes
using BlockArrays
using LinearAlgebra
# import Plots

radius = 1.0
hh=0.15
ϵ0 = 8.854e-12
μ0 = 4π*1e-7
ϵr = 2.0
μr = 3.0
c = 1/√(ϵ0*μ0)
λ = 2.9979563769321627
ω = 2π*c/λ

k0 = ω*sqrt(ϵ0*μ0)
T = boundary(CompScienceMeshes.tetmeshsphere(radius,hh))

S = BEAST.HOM(T,ϵ0*ϵr,μ0*μr,ω)
F = BEAST.FreeSpace(ϵ0,μ0,ω)
world = BEAST.create_world([2,0],[S,F])
strat = BEAST.VP()
BEAST.assign_basis!(world,strat)

lhs = BEAST.discretise_lhs(world,strat)

A(x) = [1.0,0.0,0.0]*1im/ω*exp(-1.0im*k0*x[3])
∇xA(x) = [0.0,1.0,0.0]*k0/ω*exp(-1.0im*k0*x[3])
∇A(x) = 0.0+0.0im

ex = BEAST.VPExcitation(A,∇xA,∇A,[2])
rhs = BEAST.discretise_rhs(world,ex,strat)

t =lhs[1]+lhs[2]
u = BEAST.gmres_ch(rhs[1]+rhs[2],t,world.trialdirectproductspace)

Bf = BEAST.BField(world,u[1],[1,2])
nx, nz = 100, 100
xs, zs = range(-4,stop=4,length=nx), range(-6,stop=6,length=nz)
pts = [point(0,x,z) for x in xs, z in zs]

BField = BEAST.calculate_field(pts,Bf,strat)
# display(Plots.heatmap(xs, zs, real.(getindex.(BField.field,2).+getindex.(∇xA.(pts),2)),show=true))  
# display(Plots.heatmap(xs, zs, imag.(getindex.(BField.field,2).+getindex.(∇xA.(pts),2)),show=true))  






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



X = raviartthomas(T)

κ, η = ω/c, √(μ0/ϵ0)
κ′, η′ = κ*√(ϵr*μr), η*√(μr/ϵr)

T  = Maxwell3D.singlelayer(wavenumber=κ)
T′ = Maxwell3D.singlelayer(wavenumber=κ′)
K  = Maxwell3D.doublelayer(wavenumber=κ)
K′ = Maxwell3D.doublelayer(wavenumber=κ′)

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
H = -1/(im*κ*η)*curl(E)

e = (n × E) × n
h = (n × H) × n

@hilbertspace j m
@hilbertspace k l

α, α′ = 1/η, 1/η′
pmchwt = @discretise(
    (η*T+η′*T′)[k,j] -      (K+K′)[k,m] +
         (K+K′)[l,j] + (α*T+α′*T′)[l,m] == -e[k] - h[l],
    j∈X, m∈X, k∈X, l∈X)

uu = solve(pmchwt)

import Base.Threads: @spawn
task1 = @spawn nearfield(uu[m],uu[j],X,X,κ,η,pts,E,H)
task2 = @spawn nearfield(-uu[m],-uu[j],X,X,κ′,η′,pts)

E_ex, H_ex = fetch(task1)
E_in, H_in = fetch(task2)

B_tot = H_in*μ0*μr + H_ex*μ0

# display(Plots.heatmap(xs, zs, real.(getindex.(B_tot,2))))
# display(Plots.heatmap(xs, zs, imag.(getindex.(B_tot,2))))

# maximum(norm.(B_tot-BField.field-∇xA.(pts))./norm.(B_tot))

# display(Plots.heatmap(xs,zs,norm.(B_tot-BField.field-∇xA.(pts))./maximum(norm.(B_tot))))
# display(Plots.heatmap(xs,zs,norm.(B_tot)))

# import PlotlyJS

# PlotlyJS.plot(PlotlyJS.heatmap(x=xs,y=zs,z=norm.(B_tot)))

# old = 0.008153372950474812
@test norm.(B_tot-BField.field-∇xA.(pts)))/sum(norm.(B_tot) < 0.01