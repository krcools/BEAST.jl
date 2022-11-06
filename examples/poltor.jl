using CompScienceMeshes
using LinearAlgebra

# function cone(p,q; sizemode="absolute", sizeref=2, kwargs...)
#     x = getindex.(p,1)
#     y = getindex.(p,2)
#     z = getindex.(p,3)
#     u = getindex.(q,1)
#     v = getindex.(q,2)
#     w = getindex.(q,3)
#     Plotly.cone(;x,y,z,u,v,w, sizemode, sizeref, kwargs...)
# end


fn = joinpath(dirname(pathof(CompScienceMeshes)),"geos/torus.geo")
fn = joinpath(@__DIR__, "assets/rectangular_torus.geo")
h = 0.08
Γ = CompScienceMeshes.meshgeo(fn; dim=2, h)

Γ0 = skeleton(Γ,0)
Γ1 = skeleton(Γ,1)

Σ = connectivity(Γ, Γ1, sign) 
Λ = connectivity(Γ0, Γ1, sign)

# using Plotly
# Plotly.plot([patch(Γ, opacity=0.5), wireframe(Γ)])

using BEAST
X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

K = Maxwell3D.doublelayer(gamma=0.0)
Id = BEAST.Identity()
Nx = BEAST.NCross()

M = K + 0.5Nx
qs = BEAST.defaultquadstrat(M,Y,X)
qs[2] = BEAST.DoubleNumWiltonSauterQStrat(6, 7, 6, 7, 9, 9, 9, 9)

Myx = assemble(M, Y, X, quadstrat=qs)

using LinearAlgebra
(;U,S,V) = svd(Myx)
v0 = V[:,end]

fcr, geo = facecurrents(v0, X)
# pts = [cartesian(center(chart(Γ,p))) for p in Γ]

import Plotly
pt1 = patch(Γ,norm.(fcr), opacity=0.5)
pt2 = CompScienceMeshes.cones(Γ, fcr, sizeref=0.4)
Plotly.plot([pt1, pt2])

G = assemble(Id,X,X)


norm(v0)
norm(Σ'*v0)
norm(Λ'*G*v0)

# Dict{Any, Any} with 3 entries:
#   0.5   => 0.370836
#   0.25  => 0.315732
#   0.125 => 0.265276