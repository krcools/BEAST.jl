using CompScienceMeshes, BEAST
using Base.Test

fn = joinpath(dirname(@__FILE__),"assets","sphere35.in")
m = readmesh(fn)

X = raviartthomas(m)
X1 = subset(X,1:1)

numfunctions(X)
numfunctions(X1)

t = Maxwell3D.singlelayer(gamma=im*1.0)
T1 = assemble(t,X1,X)
T2 = BEAST.assemblerow(t,X1,X)

@test T1 == T2
