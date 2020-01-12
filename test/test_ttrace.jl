using CompScienceMeshes
using BEAST
using LinearAlgebra
using Test

o, x, y, z = euclidianbasis(3)

p1 = 2x
p2 = y
p3 = 3z
p4 = o
tet = simplex(p1, p2, p3, p4)

q = abs(CompScienceMeshes.relorientation([1,2,3],[1,2,3,4]))
tri = simplex(p1, p2, p3)
n = normal(tri)

i = 3
j = 3
edg = simplex(p1, p2)

T = Float64
x = BEAST.NDLCDRefSpace{T}()
y = BEAST.NDRefSpace{T}()

p = neighborhood(tet, [0.5, 0.5, 0.0])
r = neighborhood(tri, [0.5, 0.5])

@test carttobary(edg, cartesian(p)) ≈ [0.5]
@test carttobary(edg, cartesian(r)) ≈ [0.5]

xp = x(p)[j].value
yr = y(r)[i].value

a, b = extrema((n × xp) ./ yr)
@test a ≈ b

tgt = p2 - p1
a2 = dot(n × xp, tgt) / dot(yr, tgt)
@test a ≈ a2

z = a * yr
@test z ≈ n × xp

volume(simplex(p1,p2,p4))
volume(simplex(p1,p2))

Q = BEAST.ttrace(x, tet, q, tri)

p = neighborhood(tet, [1/3, 1/3, 1/3])
r = neighborhood(tri, [1/3, 1/3])
@test cartesian(p) ≈ cartesian(r)

xp = n × x(p)[j].value
yr = y(r)[i].value
@test xp ≈ yr * Q[i,j]


x = BEAST.NDLCCRefSpace{Float64}()
y = BEAST.RTRefSpace{Float64}()
q = 3
fc = BEAST.faces(tet)[q]
Q = BEAST.ttrace(x, tet, q, fc)

nbdi = center(fc)
nbdj = neighborhood(tet, carttobary(tet, cartesian(nbdi)))

xvals = x(nbdj)
yvals = y(nbdi)

for j in 1:6
    trc = sum(Q[i,j]*yvals[i].value for i in 1:3)
    @test isapprox(trc, normal(fc) × xvals[j].value, atol=1e-4)
end
