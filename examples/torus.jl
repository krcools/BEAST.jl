using CompScienceMeshes, BEAST, BogaertInts10

fn = joinpath(dirname(@__FILE__),"assets","torus.msh")
m = CompScienceMeshes.read_gmsh_mesh(fn)

X = raviartthomas(m)
Y = buffachristiansen(m)

κ = 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
k = Maxwell3D.doublelayer(wavenumber=κ)
n = BEAST.NCross()

verts = skeleton(m,0)
edges = skeleton(m,1)
faces = skeleton(m,2)
D = connectivity(edges, faces)

M = assemble(k+0.5n,Y,X)
h = nullspace(M,0.051)
#@assert size(h,2) == 1

fcr, geo = facecurrents(h[:,end],X)

using MATLAB
mat"""
patch('Vertices',$(vertexarray(m)), 'Faces',$(cellarray(m)), 'FaceColor','flat', 'FaceVertexCData',$(norm.(fcr)))
"""
