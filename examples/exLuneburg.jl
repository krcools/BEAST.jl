using BEAST
using CompScienceMeshes
using StaticArrays
using EMwavesBEP
using Plots
using DelimitedFiles
using LinearAlgebra
using IterativeSolvers

Ω = CompScienceMeshes.tetmeshsphere(1,0.15)

Γ = boundary(Ω)
X2 = BEAST.raviartthomas(Γ)
X3 = BEAST.buffachristiansen(Γ)
X = BEAST.nedelecc3d(Ω)
ttrX = BEAST.ttrace(X,Γ)
Y = curl(X)
ttrY = BEAST.ttrace(Y,Γ)

ϵo = 1.0
μo = 1.0
ω = 7.0
κo = ω*sqrt(μo*ϵo)
pfc = im*μo*ω

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κo)
e = (n × E) × n
H = -1/(pfc)*curl(E)
h = (n × H) × n

Id = BEAST.Identity()
N = NCross()
T = Maxwell3D.singlelayer(wavenumber=κo)
K = Maxwell3D.doublelayer(wavenumber=κo)

function ω2ϵ(p)
    ω^2*(2 - (BEAST.norm(p)^2))
end

function μ1(p)
    1 #PMCHWT
end

MP = BEAST.Multiplicative(ω2ϵ)
MP2 = BEAST.Multiplicative(μ1)

A1 = assemble(MP2,Y,Y)
A2 = assemble(MP,X,X)
A3a = assemble(T,ttrX,ttrX)
A3b1 = assemble((K-0.5N),ttrX,X2)
A3b2 = assemble(T,X2,X2)
A3b3 = assemble((K+0.5N),X2,ttrX)

S = (sqrt(ϵo/μo)*A3a+sqrt(ϵo/μo)*(A3b1*inv(A3b2)*A3b3))
Aj = A1-A2-pfc*S

b1 = assemble(e,X3)
b2 = assemble(h,ttrX)

Si1 = assemble(T,ttrX,X2)*inv(Array(assemble(N,X3,X2)))
Si2 = A3b1*inv(A3b2)*assemble((K+0.5N),X2,X2)*inv(Array(assemble(N,X3,X2)))
bj = sqrt(ϵo/μo)*(Si1+Si2)*b1-b2

uj = zeros(length(bj))*im
uj = gmres!(uj, Aj, pfc*bj)
um = -uj/pfc
#um = (-S*uj-bj)
#ttrY = ttrX
#Y = X

using Plots
print("hier")
import PlotlyJS
fcrj, geo = facecurrents(uj,ttrX)
print(extrema(norm.(fcrj)))
mini,maxi = extrema(norm.(fcrj))
display(PlotlyJS.plot(patch(geo, norm.(fcrj),(mini,maxi))))
fcrm, geo = facecurrents(um,ttrY)
print(extrema(norm.(fcrm)))
display(PlotlyJS.plot(patch(geo, norm.(fcrm),(mini,maxi))))

t = range(-2,2,length=200)
t2 = range(-2,2,length=200)
pts = [point(0,s,s2) for s in t, s2 in t2]
SLj = potential(BEAST.MWSingleLayerField3D(κo),pts,um,ttrY)
SLj = reshape(SLj,(200,200))
DLm = potential(BEAST.MWDoubleLayerField3D(κo),pts,uj,ttrX)
DLm = reshape(DLm,(200,200))
Ep = E.(pts)
Hp = H.(pts)

vals = BEAST.grideval(pts, uj, X)
vals2 = BEAST.grideval(pts, um, Y)

timesteps = 100
tijd = range(0,stop=2π,length=timesteps)
vt = getindex.(vals,1)
vt2 = getindex.(vals2,1)
Ep = E.(pts)
type = Array{Complex{Float64},1}
Ep2 = vec(getindex.(reshape(Ep-DLm+SLj,(40000,1)),1))
i = 1
for pt in pts
    global i
    if pt[2]^2+pt[3]^2 < 1
        Ep2[i] = 0
    end
    i += 1
end

ϕ = range(0,stop=2π,length=100)
display(plot(surface(t,t2,getindex.(real(vt+Ep2),1),camera=(0,90),clim=(-1.5, 1.5),aspect_ratio=:equal,size=(500,400),c=ColorGradient([:red,:black,:blue]))))
display(plot!(cos.(ϕ), sin.(ϕ), ones(length(ϕ)), color="green", linewidth = 3, leg=false))

function valt(tijd)
    (real(exp.(im*tijd)*vt))
end
function valt2(tijd)
    real(exp.(im*tijd)*vt2)
end
function Et(tijd)
    (real(exp.(im*tijd)*Ep2))
end

vv = zeros(0)
vv2 = zeros(0)
Ept =zeros(0)
for time in tijd
    global vv
    global vv2
    global Ept
    append!(vv,valt(time))
    append!(vv2,valt2(time))
    append!(Ept,Et(time))
end

anim = @animate for i ∈ 1:timesteps
    plot(surface(t,t2,vv[length(vt)*(i-1)+1:length(vt)*i],camera=(0,90),clim=(-1.5, 1.5),aspect_ratio=:equal,size=(500,400),c=ColorGradient([:red,:black,:blue])))
    plot!(cos.(ϕ), sin.(ϕ), ones(length(ϕ)), color="green", linewidth = 3, leg=false)
end
gif(anim, "Luneburg_c.gif", fps = 15)

anim = @animate for i ∈ 1:timesteps
    plot(surface(t,t2,vv[length(vt)*(i-1)+1:length(vt)*i],camera=(0,90),clim=(-1.5, 1.5),aspect_ratio=:equal,size=(500,400),c=ColorGradient([:red,:black,:blue])))
end
gif(anim, "Luneburg_b.gif", fps = 15)

anim = @animate for i ∈ 1:timesteps
    plot(surface(t,t2,vv[length(vt)*(i-1)+1:length(vt)*i]+Ept[length(vt)*(i-1)+1:length(vt)*i],camera=(0,90),clim=(-1.5, 1.5),aspect_ratio=:equal,size=(500,400),c=ColorGradient([:red,:black,:blue])))
    plot!(cos.(ϕ), sin.(ϕ), ones(length(ϕ)), color="green", linewidth = 3, leg=false)
end
gif(anim, "Luneburg_a.gif", fps = 15)
