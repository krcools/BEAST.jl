using BEAST
using CompScienceMeshes
using Test

x = Vec(1.0,0.0,0.0)
y = Vec(0.0,1.0,0.0)
z = Vec(0.0,0.0,1.0)

g = BEAST.creategaussian(1.0, 3.0)
pw = BEAST.planewave(x, z, g)

t = linspace(0.0, 6.0, 100);
y = [pw([0.0,0.0,0.0],x)[1] for x in t];

Δx = 1.0
Γ = meshrectangle(10.0, Δx, Δx)
X = raviartthomas(Γ)

Δt = 1.0
Nt = 20
U = timebasisc0d1(Float64, Δt, Nt)

b = zeros(Float64, numfunctions(X), numfunctions(U))
store(v, m, k) = (b[m,k] += v)
BEAST.assemble!(pw, X, U, store)
