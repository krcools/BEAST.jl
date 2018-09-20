# Conclusion: Sauter-Schwab quadrature does not work when the trial and
# test meshes are not geometrically conforming.

using CompScienceMeshes, BEAST

fn = joinpath(dirname(@__FILE__),"./assets","torus264.msh")
m = CompScienceMeshes.read_gmsh_mesh(fn)
V, F = vertexarray(m), cellarray(m)

X = raviartthomas(m)
Y = buffachristiansen(m)

κ = 0.0
k = Maxwell3D.doublelayer(wavenumber=κ)
b = BEAST.NCross()

verts = skeleton(m,0)
edges = skeleton(m,1)
faces = skeleton(m,2)
Λ = connectivity(verts, edges)
Σᵀ = connectivity(edges, faces)
@assert all(sum(Σᵀ,dims=1) .== 0)

op = k

BEAST.quadrule(op::BEAST.MaxwellOperator3D, g::BEAST.RTRefSpace, f::BEAST.RTRefSpace, i, τ, j, σ, qd) = BEAST.qrdf(op, g, f, i, τ, j, σ, qd)
M1 = assemble(k,Y,X)
BEAST.quadrule(op::BEAST.MaxwellOperator3D, g::BEAST.RTRefSpace, f::BEAST.RTRefSpace, i, τ, j, σ, qd) = BEAST.qrss(op, g, f, i, τ, j, σ, qd)
M2 = assemble(k,Y,X)

using LinearAlgebra

G = assemble(b,Y,X)
@show norm(Σᵀ*G*Λ)

@show norm(Λ'*M1*Λ)
@show norm(Λ'*M2*Λ)

@show norm(Σᵀ*M1*Λ)
@show norm(Σᵀ*M2*Λ)

@show norm(M1)
