using CompScienceMeshes
using BEAST

fn = joinpath(dirname(@__FILE__),"assets","torus.msh")
m = CompScienceMeshes.read_gmsh_mesh(fn)

X = raviartthomas(m)
Y = buffachristiansen(m)

N = BEAST.NCross()
K = Maxwell3D.doublelayer(gamma=0.0)

verts = skeleton(m,0)
edges = skeleton(m,1)
faces = skeleton(m,2)
Λ = connectivity(verts, edges)
Σ = connectivity(faces, edges)

M = assemble(K+0.5N,Y,X)

using LinearAlgebra
s = svdvals(M)

using Plots
plot(log10.(s))

x = nullspace(M, 1.1*s[end])

norm(x)
norm(M*x)
norm(Σ'*Λ)
