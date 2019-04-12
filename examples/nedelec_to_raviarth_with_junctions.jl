using BEAST
using CompScienceMeshes

fn = joinpath(@__DIR__,"assets","thick_cladding.msh")
G01 = CompScienceMeshes.read_gmsh_mesh(fn, physical="Gamma01")
G02 = CompScienceMeshes.read_gmsh_mesh(fn, physical="Gamma02")
G12 = CompScienceMeshes.read_gmsh_mesh(fn, physical="Gamma12")

G = CompScienceMeshes.read_gmsh_mesh(fn)

@assert numcells(G01) + numcells(G02) + numcells(G12) == numcells(G)

# Build region boundaries with normals pointing inward:
G1 = weld(-G01, -G12)
@assert CompScienceMeshes.isoriented(G1)

G2 = weld(-G02, G12)
@assert CompScienceMeshes.isoriented(G2)

E1 = CompScienceMeshes.interior(G1)
E2 = CompScienceMeshes.interior(G2)
E = CompScienceMeshes.interior(G)

Nd = BEAST.nedelec(G, E)

Nd1 = BEAST.nedelec(G1, E1)
Nd2 = BEAST.nedelec(G2, E2)

RT1 = n × Nd1
RT2 = n × Nd2

A1 = CompScienceMeshes.embedding(E1,E)
A2 = CompScienceMeshes.embedding(E2,E)
