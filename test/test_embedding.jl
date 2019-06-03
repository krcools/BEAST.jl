using Test
using CompScienceMeshes
using BEAST

# G0 is the interior of the metal (no simulation required here)
# G1 is the unbounded exterior
# G2 is the bounded interior

fn = joinpath(@__DIR__, "assets", "thick_cladding.msh")
G01 = CompScienceMeshes.read_gmsh_mesh(fn, physical="Gamma01")
G02 = CompScienceMeshes.read_gmsh_mesh(fn, physical="Gamma02")
G12 = CompScienceMeshes.read_gmsh_mesh(fn, physical="Gamma12")

G = CompScienceMeshes.read_gmsh_mesh(fn)

@test numcells(G01) + numcells(G02) + numcells(G12) == numcells(G)

# Build region boundaries with normals pointing inward:
G1 = weld(-G01, -G12)
@test CompScienceMeshes.isoriented(G1)

G2 = weld(-G02, G12)
@test CompScienceMeshes.isoriented(G2)

isclosed(m) = (numcells(boundary(m)) == 0)

@test isclosed(G1)
@test isclosed(G2)

# Test the validity of the embedding matrices
E = skeleton(G,1)
E1 = skeleton(G1,1)
A1 = CompScienceMeshes.embedding(E1,E)

nd = BEAST.nedelec(G,E)
@test numfunctions(nd) == numcells(E)
nd1 = BEAST.nedelec(G1,E1)
@test numfunctions(nd1) == numcells(E1)
Z_direct = assemble(BEAST.Identity(), nd1, nd)
Z_expand = assemble(BEAST.Identity(), nd1, nd1) * A1
@test norm(Z_direct - Z_expand, Inf) ≤ √(eps(1.0))
