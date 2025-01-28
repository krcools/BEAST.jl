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

E = Maxwell3D.planewave(
    direction=point(0,0,1),
    polarization=point(1,0,0),
    wavenumber=1.0)

e = (n × E) × n

struct SingField{T} <: BEAST.Functional
    a::T
    b::T
end

function (f::SingField)(p)
    x = cartesian(p)
    x[3] ≈ 1.0 || return point(0,0,0)
    return sqrt(f.b - x[2]) * sqrt(x[2] - f.a) / sqrt(f.b-x[1]) / sqrt(x[1]-f.a) * point(0,1,0)
end

function BEAST.integrand(f::SingField, gx, ϕx)
    return dot(gx[1], ϕx)
end

e = SingField(0.0, 1.0)

b = assemble(e, X)
G = assemble(id, X, X)
u = G \ b

import Plotly
fcr, geo = facecurrents(u, X)
Plotly.plot(
    [
        # patch(geo, norm.(fcr), opacity=0.5),
        CompScienceMeshes.cones(geo, real.(fcr), sizeref=10.0),
    ])

q = adjoint(svdA.U) * HY * b
plot(log10.(abs.(q)))
plot!(log10.(abs.(svdA.S)))