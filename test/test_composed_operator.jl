######## test multiplied basis
using CompScienceMeshes
using LinearAlgebra
using BEAST
using Test
using StaticArrays


# different options to compute selfpatch trace
Γ = meshcuboid(1.0,1.0,1.0,0.4; generator=:compsciencemeshes)
Γd = BEAST.GlobalDisplacementMesh(Γ,1.0)
X = raviartthomas(Γ)
Xd = raviartthomas(Γd)

∇G = BEAST.HH3DGradGreen(0.0)
Z = BEAST.PotentialIntegralOperator{2}(∇G,×,b->b)
K = BEAST.ttrace(Z,false;testfunction_tangential=false) #-1.0 trace taken allong the test normal (from outside to inside)
M = assemble(K,X,X;quadstrat = [BEAST.SingleNumQStrat(6),BEAST.DoubleNumSauterQstrat(6,6,6,6,6,6)])



K_compare = Maxwell3D.doublelayer(wavenumber = 0.0) - 0.5*NCross() #trace from outside to inside
M_compare = assemble(K_compare,X,X;quadstrat = [BEAST.SingleNumQStrat(6),BEAST.DoubleNumSauterQstrat(6,6,6,6,6,6)])

M_displaced = assemble(K,Xd,X;quadstrat = [BEAST.SingleNumQStrat(6),BEAST.DoubleNumSauterQstrat(6,6,6,6,6,6)])

@test M_displaced ≈ M
@test M_displaced ≈ M_compare 

####################################################
#
# Material interactions (global multi-trace context)
#
####################################################

Γ_mirror = -mirrormesh(Γ,(@SVector [1.0,0.0,0.0]),(@SVector [0.0,0.0,0.0]))
Γd_mirror = BEAST.GlobalDisplacementMesh(Γ_mirror,-1.0)
Xm = raviartthomas(Γ_mirror)
Xdm = raviartthomas(Γd_mirror)

K = ttrace(Z,true) #number does not matter here, trace part should just be spawned, this number only matters fr selfpatch or if user is not carefull with displacement meshes
M = assemble(K,Xm,X;quadstrat = [BEAST.SingleNumQStrat(6),BEAST.DoubleNumSauterQstrat(6,6,6,6,6,6)])

K_compare = Maxwell3D.doublelayer(wavenumber = 0.0) - 0.5*NCross() #trace from outside to inside
M_compare = assemble(K_compare,Xm,X;quadstrat = [BEAST.SingleNumQStrat(6),BEAST.DoubleNumSauterQstrat(6,6,6,6,6,6)])

M_displaced = assemble(K,Xdm,X;quadstrat = [BEAST.SingleNumQStrat(6),BEAST.DoubleNumSauterQstrat(6,6,6,6,6,6)])

@test M_displaced ≈ M
@test M_displaced ≈ M_compare 



