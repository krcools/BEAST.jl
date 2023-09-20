using BEAST
using StaticArrays
using CompScienceMeshes
using Test
k = 1.0

green = HHH.green(wavenumber=k)
gradgreen = HHH.gradgreen(wavenumber=k)
b = basisfunction()
green
gradgreen
n×green
n×gradgreen
gradgreen×nothing
gradgreen×n
gradgreen(n×b)
gradgreen(n*b)
n×gradgreen(n×b)


Γ1 = boundary(Mesh([SVector{3}([0.0,0.0,0.0]),SVector{3}([1.0,0.0,0.0]),SVector{3}([0.0,1.0,0.0]),SVector{3}([0.0,0.0,1.0])],[SVector{4}([1,2,3,4])]))
Γ2 = CompScienceMeshes.translate(Γ1,[0.0,0.0,0.0])

X1 = raviartthomas(Γ1)
X2 = raviartthomas(Γ2)
X3 = lagrangec0d1(Γ1) 

m1 = assemble(n×(n×(gradgreen×nothing)(n×(n×b))),X1,X2)
m2 = assemble((gradgreen)(n⋅(n×b)),X1,X2)
m3 = assemble(n⋅(gradgreen)(n⋅(n×b)),X3,X2)
m4 = assemble((green)(n×(n*b)),X1,X3)
dl = assemble(Maxwell3D.doublelayer(wavenumber=k),X1,X2)

@test abs(sum(m1-dl)+sum(m2)+sum(m3)+sum(m4)) < 10^(-16)