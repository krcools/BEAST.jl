using BEAST, CompScienceMeshes

#=
example file for which the nonconforming assembly
strategy fails
=#

mesh_h = 0.2/5
h_ref = 0.2/10
m = meshcuboid(1.0,1.0,0.25,mesh_h)
mref = meshcuboid(1.0,1.0,0.25,h_ref)

Xref = BEAST.raviartthomas(mref)
X = BEAST.raviartthomas(m)
t = Maxwell3D.singlelayer(wavenumber=1.0)

default_qstrat = BEAST.defaultquadstrat(t, X, Xref)
qstrat = BEAST.NonConformingIntegralOpQStrat(default_qstrat)
Z = BEAST.assemble(t, X, Xref, quadstrat=qstrat)

#element 2725 in m and element 10673 in mref are problematic pairs
import Plotly
m1 = CompScienceMeshes.Mesh(m.vertices, [m.faces[2725]])
m2 = CompScienceMeshes.Mesh(mref.vertices, [mref.faces[10673]])
Plotly.plot([wireframe(m1), wireframe(m2)])

s1 = CompScienceMeshes.simplex(m1.vertices[m1.faces[1]])
s2 = CompScienceMeshes.simplex(m2.vertices[m2.faces[1]])
CompScienceMeshes.overlap(s1, s2)

