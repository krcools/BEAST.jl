##### user manual composed operator and potential

using BEAST
using CompScienceMeshes

Γ = meshcuboid(1.0,1.0,1.0,0.2);
X = lagrangecxd0(Γ)

S = BEAST.PotentialIntegralOperator{2}(BEAST.HH3DGreen(1im),*,B->B)
St = BEAST.trace(S,1.0)
assemble(St,X,X)

SB = Helmholtz3D.singlelayer(wavenumber=1.0)
assemble(SB,X,X)

Y = raviartthomas(Γ)
K = BEAST.PotentialIntegralOperator{2}(BEAST.HH3DGradGreen(0*1im),×,B->B)
Kt = BEAST.ttrace(K,1.0;testfunction_tangential=true)
assemble(Kt,Y,Y)

KB = Maxwell3D.doublelayer(wavenumber=0.0)