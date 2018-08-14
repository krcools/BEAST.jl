using CompScienceMeshes
using BEAST
using Test

p1 = point(0,0,0)
p2 = point(1/3,0,0)
p3 = point(0,1/3,0)

q1 = point(0,0,0)
q2 = point(0.125,0,0)

m = Mesh([p1,p2,p3], [index(1,2,3)])
n = Mesh([q1,q2], [index(1,2)])
translate!(n, point(0,0,20))

X = lagrangec0d1(m, boundary(m))
@test numfunctions(X) == 3

x = refspace(X)
s = chart(m, m.faces[X.fns[1][1].cellid])
c = neighborhood(s, [1,1]/3)      # get the barycenter of that patch
v = x(c, Val{:withcurl})       # evaluate the Lagrange elements in c, together with their curls

@test (s[3] - s[2]) / (2 * volume(s)) ≈ v[1][2]
@test (s[1] - s[3]) / (2 * volume(s)) ≈ v[2][2]
@test (s[2] - s[1]) / (2 * volume(s)) ≈ v[3][2]

Y = lagrangec0d1(n, boundary(n))
@test numfunctions(Y) == 2

y = refspace(Y)
t = chart(n, n.faces[Y.fns[1][1].cellid])
d = neighborhood(t, [1]/2)
tg = normalize(tangents(d,1))
w = y(d)

@test w[1][1] ≈ 0.5
@test w[2][1] ≈ 0.5

κ = 0.0
T = BEAST.NitscheHH3(κ)
Tyx = assemble(T, Y, X)

@test size(Tyx) == (numfunctions(Y), numfunctions(X))

R = norm(cartesian(c)-cartesian(d))
estimate = volume(t) * volume(s) * dot(v[1][2], w[1][1] * tg) / (4π*R)
actual = Tyx[1,1]

@test (estimate - actual) / actual < 0.005

## test the value of the skeleton gram matrix
I = BEAST.Identity()
Iyy = assemble(I, Y, Y)

@test size(Iyy) == (2,2)

qps = BEAST.quadpoints(y, [t], (10,))[1,1]
estimate = 0.0
for qp in qps
    igd = qp.value[1][1] * qp.value[2][1]
    global estimate += qp.weight * igd
end
actual = Iyy[1,2]
@test norm(estimate - actual) < 1e-6
