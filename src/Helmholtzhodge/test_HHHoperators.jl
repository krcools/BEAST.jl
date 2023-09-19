using BEAST
using StaticArrays
using CompScienceMeshes
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

display(assemble(gradgreen×nothing,X1,X2))
display(assemble(Maxwell3D.doublelayer(wavenumber=k),X1,X2))