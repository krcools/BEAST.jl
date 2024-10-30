using Test

using CompScienceMeshes
using BEAST

s = simplex(
    point(5,1,0),
    point(-1,3,0),
    point(0,0,2),
)

locspace = BEAST.RTRefSpace{Float64}()
t = Maxwell3D.singlelayer(wavenumber=1.0)

igd = BEAST.Integrand(t, locspace, locspace, s, s)

u = point(1/3,1/3)
v = point(1/4,2/3)

igd(u,v)

dom = CompScienceMeshes.domain(s)
p = CompScienceMeshes.neighborhood(dom, u)
f̂ = locspace(p)
Dx = tangents(s, u)
f = map(f̂) do fi
    (value = Dx * fi.value, divergence = fi.divergence) end

q = neighborhood(s, u)
g = locspace(q)
J = jacobian(q)