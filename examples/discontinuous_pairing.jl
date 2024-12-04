# compute the continuity and infsup constants for various continuous and discontinuous pairings

using CompScienceMeshes
using BEAST

using Makeitso
using BakerStreet

h = 0.1
Γ = meshcuboid(1.0, 1.0, 1.0, h)

id = BEAST.Identity()
nx = BEAST.NCross()
s = Maxwell3D.singlelayer(gamma=1.0)
t = Maxwell3D.singlelayer(wavenumber=1.0)

X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

A = assemble(id, Y, X)
B = assemble(nx, Y, X)
C = assemble(t, X, X)

SXX = Symmetric(-assemble(s, X, X))
SYY = Symmetric(-assemble(s, Y, Y))

using LinearAlgebra
HX = sqrt(SXX)
HY = sqrt(SYY)

QX = inv(HX)
QY = inv(HY)

svdA = svd(QY*A*QX)
svdB = svd(QY*B*QX)
svdC = svd(QX*C*QX)

using Plots
plot(svdA.S)
plot!(svdB.S)
plot!(svdC.S)
