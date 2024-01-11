using CompScienceMeshes, BEAST, LinearAlgebra
#Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
Γ = meshsphere(radius=1.0, h=0.30)

X = lagrangecxd0(Γ)

Δt = 0.16
Nt = 300
T = timebasisshiftedlagrange(Δt, Nt, 0)
U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U

width, delay, scaling = 16.0, 24.0, 1.0
gaussian = creategaussian(width, delay, scaling)
e = BEAST.planewave(point(0,0,1), 1.0, gaussian)

@hilbertspace j
@hilbertspace j′

SL = TDAcustic3D.acusticsinglelayer(speedofsound=1.0, numdiffs=1)
# BEAST.@defaultquadstrat (SL, W, V) BEAST.OuterNumInnerAnalyticQStrat(7)

tdacusticsl = @discretise SL[j′,j] == -1.0e[j′]   j∈V  j′∈W
xacusticsl = solve(tdacusticsl)

Z = assemble(SL,W,V)
#corregere da qui 
#import Plots
#Plots.plot(xacusticsl[1,250:300])

#import Plotly
#fcr, geo = facecurrents(xefie[:,125], X)
#Plotly.plot(patch(geo, norm.(fcr)))

import BEAST.ConvolutionOperators

for a in 0:300
    za=ConvolutionOperators.timeslice(Z,a)
    #zawilton=ConvolutionOperators.timeslice(Zs,a)
    for i in 1:numfunctions(X)
        for j in 1:numfunctions(X)
            if  (za[i,j])<0.0
                println(a,[i,j]," ",za[i,j]," ",)
            end
        end
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

