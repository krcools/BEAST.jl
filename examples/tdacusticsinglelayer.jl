using CompScienceMeshes, BEAST, LinearAlgebra
#Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
Γ = CompScienceMeshes.meshsphere(1.0,0.3)#CompScienceMeshes.meshcuboid(1.0,1.0,1.0,0.3)

X = lagrangecxd0(Γ)

Δt, Nt = 0.08, 600 
T = timebasisshiftedlagrange(Δt, Nt, 0)
U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U

width, delay, scaling = 16.0, 24.0, 1.0
gaussian = creategaussian(width, delay, scaling)
e = BEAST.planewave(point(0,0,1), 1.0, gaussian)

@hilbertspace j
@hilbertspace j′

SL = TDAcustic3D.acusticsinglelayer(speedofsound=1.0, numdiffs=0)
# BEAST.@defaultquadstrat (SL, W, V) BEAST.OuterNumInnerAnalyticQStrat(7)

tdacusticsl = @discretise SL[j′,j] == -1.0e[j′]   j∈V  j′∈W
xacusticsl = solve(tdacusticsl)

Z = assemble(SL,W,V)
#corregere da qui 
import Plots

Plots.plot(xacusticsl[1,1:300],label="Current")
Plots.xlabel!("t")
Plots.savefig("stablecurrent.png") 

xacusticsl[1,200:300]
#pval=ConvolutionOperators.polyvals(Z)
import Plotly
#fcr, geo = facecurrents(xefie[:,125], X)
#Plotly.plot(patch(geo, norm.(fcr)))

import BEAST.ConvolutionOperators

for a in 1:1
    za=ConvolutionOperators.timeslice(Z,a)
    zawilton=ConvolutionOperators.timeslice(Zs,a)
    for i in 1:10 #numfunctions(X)
        #for j in 1:numfunctions(X)
            if  norm(za[i,i]-zawilton[i,i])>=10^(-6)
                println(a," ",i," ",za[i,i]," ",zawilton[i,i])
            end
        #end
    end
end



name=readline()
parse(Float64, name)

Xefie, Δω, ω0 = fouriertransform(xefie, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω.-1.0))

ω1 = ω[i1]
ue = Xefie[:,i1] / fouriertransform(gaussian)(ω1)

fcr, geo = facecurrents(ue, X)
Plotly.plot(patch(geo, norm.(fcr)))

