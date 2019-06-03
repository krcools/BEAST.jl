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
G = weld(G,G12,seam=b)
@test numcells(G01) + numcells(G02) + numcells(G12) == numcells(G)
@test numcells(skeleton(G01,1)) + numcells(skeleton(G02,1)) + numcells(skeleton(G12,1)) == numcells(skeleton(G,1)) - 2*numcells(b)
