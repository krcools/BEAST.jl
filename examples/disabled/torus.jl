using CompScienceMeshes, BEAST

fn = joinpath(dirname(pathof(BEAST)),"../examples/assets/torus.msh")
m = CompScienceMeshes.read_gmsh_mesh(fn)
@show numcells(m)

X = raviartthomas(m)
@show numfunctions(X)
Y = buffachristiansen(m)
@show numfunctions(Y)

κ = 1.0e-5
S = Maxwell3D.singlelayer(wavenumber=κ)
D = Maxwell3D.doublelayer(wavenumber=κ)
N = BEAST.NCross()

verts = skeleton(m,0)
edges = skeleton(m,1)
faces = skeleton(m,2)
Λ = connectivity(verts, edges)
Σᵀ = connectivity(edges, faces)
@assert all(sum(Σᵀ,dims=1) .== 0)

Sxx = assemble(S,X,X)
Myx = assemble(D+0.5N,Y,X)

using LinearAlgebra

norm(Σᵀ*Myx*Λ)

ϵ, μ = 1.0, 1.0
ω = κ / √(ϵ*μ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
H = -1/(im*μ*ω)*curl(E)
e = (n × E) × n
h = (n × H) × n

ex = assemble(e, X)
hy = assemble(h, Y)

u1 = Sxx \ ex
u2 = Myx \ hy

Φ, Θ = [0.0], range(0,stop=π,length=100)
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]
ff1 = potential(MWFarField3D(wavenumber=κ), pts, u1, X)
ff2 = potential(MWFarField3D(wavenumber=κ), pts, u2, X)

# Now for the classic MFIE
C = BEAST.DoubleLayerRotatedMW3D(im*κ)
J = BEAST.Identity()

Mxx = assemble(C-0.5J,X,X)
hx = assemble(n × H,X)
u3 = Mxx \ hx
ff3 = potential(MWFarField3D(wavenumber=κ), pts, u3, X)

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
