using Test
using CompScienceMeshes, BEAST

fn = joinpath(dirname(@__FILE__),"assets","sphere35.in")
m = readmesh(fn)

τ = chart(m,cells(m)[1])
σ = chart(m,cells(m)[end])

X = lagrangecxd0(m)
x = refspace(X)

S = BEAST.Helmholtz3D.singlelayer(gamma=1.5im)
qd = BEAST.quaddata(S, x, x, [τ,σ], [τ,σ])
