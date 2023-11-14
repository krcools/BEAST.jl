using BEAST
using CompScienceMeshes
using StaticArrays
using Test

tcell = simplex((@SVector [1.0,0.0,0.0]), (@SVector [0.0,1.0,0.0]),(@SVector [0.0,0.0,0.0]))
bcell = simplex((@SVector [0.0,0.8,0.0]),(@SVector [1.0,1.0,0.0]),(@SVector [0.9,0.0,0.0]))

bshapes = BEAST.RTRefSpace{Float64}()
tshapes = BEAST.LagrangeRefSpace{Float64,1,3,3}()

biop = BEAST.DotLocal(nt,(BEAST.DoubleIntegralR{Float64}()×B))

quadstrat = BEAST.DoubleNumSauterQstrat(6,7,6,6,6,6) 

qd = quaddata(biop, tshapes, bshapes, [tcell], [bcell], quadstrat)
zlocal1 = zeros(Float64, 6, 6)
zlocal2 = zeros(Float64, 6, 6)
##without subdevision
qrule = quadrule(biop, tshapes, bshapes, 1, tcell, 1, bcell, qd, quadstrat)
BEAST.momintegrals!(biop, tshapes, bshapes, tcell, bcell, zlocal1, qrule)

##with subdevision
BEAST.momintegrals_overlap!(biop, tshapes, bshapes, tcell, bcell, zlocal2, quadstrat)
@test zlocal1 ≈ zlocal2