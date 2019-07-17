using LinearAlgebra
using Test

using CompScienceMeshes
using BEAST

import CompScienceMeshes.SComplex2D

fn = joinpath(@__DIR__, "assets", "thick_cladding.msh")
G01 = SComplex2D(CompScienceMeshes.read_gmsh_mesh(fn, physical="Gamma01"))
G02 = SComplex2D(CompScienceMeshes.read_gmsh_mesh(fn, physical="Gamma02"))
G12 = SComplex2D(CompScienceMeshes.read_gmsh_mesh(fn, physical="Gamma12"))

G = SComplex2D(CompScienceMeshes.read_gmsh_mesh(fn))

b = boundary(G02)
G = weld(G01,G02,seam=b)
@test numcells(G01) + numcells(G02) == numcells(G)
@test numcells(skeleton(G01,1)) + numcells(skeleton(G02,1)) == numcells(skeleton(G,1)) + numcells(b)

G = weld(G,G12,seam=b)
@test numcells(G01) + numcells(G02) + numcells(G12) == numcells(G)
@test numcells(skeleton(G01,1)) + numcells(skeleton(G02,1)) + numcells(skeleton(G12,1)) == numcells(skeleton(G,1)) + 2*numcells(b)

G1 = weld(-G01, -G12, seam=b)
@test CompScienceMeshes.isoriented(G1)
G2 = weld(-G02, G12, seam=b)
@test CompScienceMeshes.isoriented(G2)

isclosed(m) = (numcells(boundary(m)) == 0)

@test isclosed(G1)
@test isclosed(G2)

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
