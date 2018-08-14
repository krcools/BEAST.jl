using CompScienceMeshes
using BEAST

using Test

m1 = meshrectangle(1.0, 1.0, 0.5)
m2 = translate!(meshrectangle(1.0, 1.0, 0.5), point(0,0,1))
m3 = translate!(meshrectangle(1.0, 1.0, 0.5), point(0,0,2))

X1 = raviartthomas(m1)
X2 = raviartthomas(m2)
X3 = raviartthomas(m3)

X = X1 × X2 × X3

n1 = numfunctions(X1)
n2 = numfunctions(X2)
n3 = numfunctions(X3)

@test numfunctions(X) == n1 + n2 + n3
nt = numfunctions(X)

T = MWSingleLayer3D(1.0)
t = assemble(T, X, X)

@test size(t) == (nt,nt)
