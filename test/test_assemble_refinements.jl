using StaticArrays
using CompScienceMeshes
using BEAST
# using SauterSchwabQuadrature
using LinearAlgebra
using Test

function generate_refpair(;angle=180)
    vertices = [
        SVector(0.0, 1.0, 0.0),
        SVector(-0.5, 0.0, 0.0),
        SVector(+0.5, 0.0, 0.0),
        SVector(0.0, 1.0*cos(angle/180*π), 1.0*sin(angle/180*π))
    ]
    
    triangles = [
        SVector(1, 2, 3),
        SVector(2, 4, 3)
    ]
    return Mesh(vertices, triangles)
end
Γ = generate_refpair(angle=45)

X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

Kop = Maxwell3D.doublelayer(wavenumber=0.0)

qs = BEAST.defaultquadstrat(Kop, X, Y)
qs = BEAST.DoubleNumWiltonSauterQStrat(1, 1, 12, 13, 12, 12, 12, 12)
@show qs
@show qs(Kop, X, Y)
@show qs(Kop, Y, X)
K1 = assemble(Kop, X, Y; quadstrat=qs)
K2 = assemble(Kop, Y, X; quadstrat=qs)

@show K1[1,1]
@show K2[1,1]

K1[1,1] - K2[1,1]
@test K1[1,1] ≈ K2[1,1] rtol=1e-7

# Verify that quadrature strategy is passed through
qds(i) = BEAST.DoubleNumWiltonSauterQStrat(3,3,8,8,i,i,i,i)
Kxy_1 = (assemble(Kop, X, Y, quadstrat=qds(1)))[1,1]
Kxy_2 = (assemble(Kop, X, Y, quadstrat=qds(5)))[1,1]
Kxy_3 = (assemble(Kop, X, Y, quadstrat=qds(10)))[1,1]
# Kxy_4 = (assemble(Kop, X, Y, quadstrat=qds(15)))[1,1]

@test 0.1 < abs(Kxy_1 - Kxy_3)/abs(Kxy_3) < 0.4
@test 0.0008 < abs(Kxy_2 - Kxy_3)/abs(Kxy_3) < 0.002

Kyx_1 = (assemble(Kop, Y, X, quadstrat=qds(1)))[1,1]
Kyx_2 = (assemble(Kop, Y, X, quadstrat=qds(5)))[1,1]
Kyx_3 = (assemble(Kop, Y, X, quadstrat=qds(10)))[1,1]
# Kyx_4 = (assemble(Kop, Y, X, quadstrat=qds(15)))[1,1]

@test 0.2 < abs(Kyx_1 - Kyx_3)/abs(Kyx_3) < 0.3
@test 0.0006 < abs(Kyx_2 - Kyx_3)/abs(Kyx_3) < 0.0007

Γtr = CompScienceMeshes.translate(Γ, SVector(2.0, 0.0, 0.0))

Γ2 = weld(Γ, Γtr)

X2 = raviartthomas(Γ2)
Y2 = buffachristiansen(Γ2)

qds2(i) = BEAST.DoubleNumWiltonSauterQStrat(i,i,i,i,10,10,10,10)
K2xy_1 = (assemble(Kop, X2, Y2, quadstrat=qds2(1)))[1,2]
K2xy_2 = (assemble(Kop, X2, Y2, quadstrat=qds2(5)))[1,2]
K2xy_3 = (assemble(Kop, X2, Y2, quadstrat=qds2(10)))[1,2]
@test 0.06 < abs(K2xy_1 - K2xy_3)/abs(K2xy_3) < 0.07
@test 1e-6 < abs(K2xy_2 - K2xy_3)/abs(K2xy_3) < 1e-5

K2yx_1 = (assemble(Kop, Y2, X2, quadstrat=qds2(1)))[1,2]
K2yx_2 = (assemble(Kop, Y2, X2, quadstrat=qds2(5)))[1,2]
K2yx_3 = (assemble(Kop, Y2, X2, quadstrat=qds2(10)))[1,2]
@test 0.06 < abs(K2yx_1 - K2yx_3)/abs(K2yx_3) < 0.07
@test 1e-6 < abs(K2yx_2 - K2yx_3)/abs(K2yx_3) < 1e-5