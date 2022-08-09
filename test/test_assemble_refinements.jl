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

K1 = assemble(Kop, X, Y)
K2 = assemble(Kop, Y, X)

@test K1[1,1] ≈ K2[1,1] rtol=1e-3