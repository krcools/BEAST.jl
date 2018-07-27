using CompScienceMeshes, BEAST, BogaertInts10

# function BEAST.quadrule(op::BEAST.MaxwellOperator3D, g::BEAST.RTRefSpace, f::BEAST.RTRefSpace, i, τ, j, σ, qd)
#
#     dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))
#     xtol = 0.2
#     k = norm(op.gamma)
#     hits = 0
#     xmin = xtol
#     for t in τ.vertices
#         for s in σ.vertices
#             d = norm(t-s)
#             xmin = min(xmin, k*d)
#             if d < dtol
#                 hits +=1
#                 break
#             end
#         end
#     end
#
#     hits == 3   && return BogaertInts10.BogaertSelfPatchStrategy(5)
#     hits == 2   && return BogaertInts10.BogaertEdgePatchStrategy(8, 4)
#     hits == 1   && return BogaertInts10.BogaertPointPatchStrategy(2, 3)
#     xmin < xtol && return BEAST.WiltonSEStrategy(
#         qd.tpoints[1,i],
#         BEAST.DoubleQuadStrategy(
#             qd.tpoints[2,i],
#             qd.bpoints[2,j], #         ),
#     )
#     return BEAST.DoubleQuadStrategy(
#         qd.tpoints[1,i],
#         qd.bpoints[1,j],
#     )
# end


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
