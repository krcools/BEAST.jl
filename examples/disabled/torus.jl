using CompScienceMeshes, BEAST

fn = joinpath(dirname(@__FILE__),"../assets","torus.msh")
m = CompScienceMeshes.read_gmsh_mesh(fn)

X = raviartthomas(m)
Y = buffachristiansen(m)

κ = 0.0
t = Maxwell3D.singlelayer(wavenumber=κ)
k = Maxwell3D.doublelayer(wavenumber=κ)
b = BEAST.NCross()

verts = skeleton(m,0)
edges = skeleton(m,1)
faces = skeleton(m,2)
Λ = connectivity(verts, edges)
Σᵀ = connectivity(edges, faces)
all(sum(Σᵀ,dims=1) .== 0)

op = k

M = assemble(k+0.5b,Y,X)

using LinearAlgebra

norm(Σᵀ*M*Λ)

G = assemble(b, Y, X)
norm(M)
norm(Σᵀ*G*Λ)


# s = svdvals(M)
# h = nullspace(M,s[end]*1.0001)
#
# fcr, geo = facecurrents(h[:,end],X)
#
# V,F = vertexarray(m), cellarray(m)
# B = [real(f[i]) for f in fcr, i in 1:3]

# using MATLAB
# mat"mesh = Mesh($V,$F)"
# mat"[C,N] = faceNormals(mesh)"
# mat"figure; hold on"
# mat"patch(mesh,$(norm.(fcr)))"
# mat"quiver3x(C,$B)"
