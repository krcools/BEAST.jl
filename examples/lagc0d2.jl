using CompScienceMeshes
using BEAST
using LinearAlgebra

import Plotly

m = meshsphere(radius=1.0, h=0.25)
nodes = skeleton(m,0)
edges = skeleton(m,1)

X = BEAST.lagrangec0d2(m, nodes, edges)
Y = BEAST.lagrangecxd0(m)

uⁱ = Helmholtz3D.planewave(wavenumber=1.0, direction=point(0,0,1))
f = strace(uⁱ,m)
fy = assemble(f,Y)

Id = BEAST.Identity()
qs = BEAST.SingleNumQStrat(8)
Gxy = assemble(Id, X, Y, quadstrat=qs)
Gxx = assemble(Id, X, X, quadstrat=qs)
Gyy = assemble(Id, Y, Y, quadstrat=qs)

Pxy = inv(Matrix(Gxx)) * Gxy * inv(Matrix(Gyy))

fX = Pxy * fy
fcr, geo = facecurrents(fX, X)
Plotly.plot(patch(m, real.(fcr)))