using CompScienceMeshes
using BEAST
using Test
using LinearAlgebra

@warn "CompScienceMeshes being patched..."
import StaticArrays
function CompScienceMeshes.chart(Smesh::CompScienceMeshes.subdMesh,E)
    element = Smesh.elements[E]
    Svertices = Smesh.vertices
    nodes = element.RingNodes
    N = length(nodes)
    verticecoords = Vector{StaticArrays.SVector{3,Float64}}(undef,N)

    for i = 1:N
        verticecoords[i] = Svertices[nodes[i]].Coords
    end
    return chart = CompScienceMeshes.subd_chart(E,N,nodes,verticecoords)
end

for T in [Float32, Float64]
G = readmesh(joinpath(dirname(@__FILE__),"assets","sphere872.in"),T=T)
 # G0 = readmesh(fn)
# G0 = meshsphere(1.0, 0.5)
#G = readmesh("/Users/Benjamin/Documents/sphere.in")
# G = Loop_subdivision(G0)
local X = subdsurface(G)

#els, ad = BEAST.assemblydata(X)

identityop  = Identity()
singlelayer = Helmholtz3D.singlelayer(gamma=T(1.0))
I = assemble(identityop, X, X)
#S = assemble(singlelayer, X, X)
ncd = cond(Matrix(I))

@test ncd ≈ T(64.50401358713235) rtol=sqrt(eps(T))
end