using CompScienceMeshes
using LinearAlgebra

function cone(p,q; sizemode="absolute", sizeref=2, kwargs...)
    x = getindex.(p,1)
    y = getindex.(p,2)
    z = getindex.(p,3)
    u = getindex.(q,1)
    v = getindex.(q,2)
    w = getindex.(q,3)
    Plotly.cone(;x,y,z,u,v,w, sizemode, sizeref, kwargs...)
end


fn = joinpath(dirname(pathof(CompScienceMeshes)),"geos/torus.geo")
Γ = CompScienceMeshes.meshgeo(fn; dim=2, h=0.5)

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
pts = [cartesian(center(chart(Γ,p))) for p in Γ]

import Plotly
pt1 = patch(Γ,norm.(fcr), opacity=0.5)
pt2 = cone(pts, fcr, sizeref=0.2)
Plotly.plot([pt1, pt2])

G = assemble(Id,X,X)


norm(v0)
norm(Σ'*v0)
norm(Λ'*G*v0)

# Study the kernel of the PMCHWT
SL = Maxwell3D.singlelayer(wavenumber=0.01)
DL = Maxwell3D.doublelayer(wavenumber=0.01)

@hilbertspace m j
@hilbertspace k l
a =
    DL[k,m] - SL[k,j] +
    SL[l,m] + DL[l,j]

A = assemble(@discretise(a, m∈X, j∈X, k∈X, l∈X))
M = Matrix(A)
(;U,S,V) = svd(M)

using Plots
plotly()
plot()

v0 = V[:,end]
v1 = V[:,end-1]
