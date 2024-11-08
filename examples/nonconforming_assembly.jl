using BEAST, CompScienceMeshes
using Infiltrator
using StaticArrays

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

default_qstrat = BEAST.defaultquadstrat(t, Xref, X)
qstrat = BEAST.NonConformingIntegralOpQStrat(default_qstrat)
Z = BEAST.assemble(t, Xref, X, quadstrat=qstrat);
error()
#element 2725 in m and element 10673 in mref are problematic pairs
# import Plotly
# m1 = CompScienceMeshes.Mesh(m.vertices, [m.faces[2725]])
# m2 = CompScienceMeshes.Mesh(mref.vertices, [mref.faces[10673]])
# Plotly.plot([wireframe(m1), wireframe(m2)])

# s1 = CompScienceMeshes.simplex(m1.vertices[m1.faces[1]])
# s2 = CompScienceMeshes.simplex(m2.vertices[m2.faces[1]])
# CompScienceMeshes.overlap(s1, s2)

tchart = simplex(
    point(0.8608841985496505, 0.0, 0.1467009531739533),
    point(0.8609389984288958, 0.0, 0.1467010505124744),
    point(0.8557537733723628, 0.0, 0.15543618661244227))

bchart = simplex(
    point(0.8605621731948716, 0.0, 0.1472492443096586),
    point(0.850767750963251, 0.0, 0.16383574186148747),
    point(0.8412930124846245, 0.0, 0.14707668556921702))
    
@show volume(tchart), volume(bchart)
@show BEAST._numhits(tchart, bchart)

# function CompScienceMeshes.Mesh(s::CompScienceMeshes.Simplex)
#     vertices = s.vertices
#     faces = [SVector(1,2,3)]
#     return CompScienceMeshes.Mesh(vertices, faces)
# end

# function CompScienceMeshes.patch(s::CompScienceMeshes.Simplex; kwargs...)
#     return CompScienceMeshes.patch(Mesh(s); kwargs...)
# end

# p1 = patch(tchart; color=:red, opacity=1.0)
# p2 = patch(bchart; color=:blue, opacity=0.25)
# Plotly.plot([p1,p2])

Xref_local = BEAST.refspace(Xref)
X_local = BEAST.refspace(X)

qd = BEAST.quaddata(t, Xref_local, X_local, [tchart], [bchart], qstrat)
qr = BEAST.quadrule(t, Xref_local, X_local, 1, tchart, 1, bchart, qd, qstrat)

# error()


out = zeros(ComplexF64,3,3)
tsing = BEAST.singularpart(t)
# BEAST.momintegrals!(t, Xref_local, X_local, tchart, bchart, out, qr)

dx, x = first(qr.outer_quad_points)
# BEAST.innerintegrals!(tsing, x, Xref_local, Xref, tchart, bchart, out, qr, dx)


using WiltonInts84
cx = cartesian(x)
b1, b2, b3 = bchart.verticesinc
WiltonInts84.wiltonints(b1, b2, b3, cx, Val{1})

a = -0.01927176760943159
b = 0-1.834267662268588e-6
p = 2.7959608455489207e-6
h = 0.0
m = SVector{3}(0.008954817714265238, 0.0, -0.9999599048160402)
UB = Val{1}

WiltonInts84.segintsg(a, b, p, h, m, UB)

T = typeof(a)
z = zero(T)
ϵ = eps(T) * 1_000

h2, d = h^2, abs(h)
q2 = p^2+h2
ra, rb = sqrt(a^2+q2), sqrt(b^2+q2)

sgn = d < ϵ ? zero(T) : sign(h)
I1 = abs(p) < ϵ ? z : sgn*(atan((p*b)/(q2+d*rb)) - atan((p*a)/(q2 + d*ra)))
#j = (q2 < ϵ^2) ? (b > 0 ? log(b/a) : log(a/b)) : log(b + rb) - log(a + ra)
a2 = a^2
b2 = b^2
# if b < 0 && q2 < max(a2,b2) * (0.5e-3)^2

# both negative, both small
j1 = log(a/b) + log((1-(q2/b2)/4) / (1-(q2/a2)/4))

# generic case
j2 = log(b + rb) - log(a + ra)

# q/b small, a/b large
j3 = log((b + rb) / (a * (-0.5*q2/a2 + 0.125*(q2/a2)^2)))
# end
abs(j2-j3)
abs(j1-j3)