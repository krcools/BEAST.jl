using CompScienceMeshes
using BEAST
using Test

p1 = point(0,0,0)
p2 = point(1/3,0,0)
p3 = point(0,1/3,0)

q1 = point(0,0,0)
q2 = point(0.125,0,0)

m = Mesh([p1,p2,p3], [CompScienceMeshes.SimplexGraph(1,2,3)])
n = Mesh([q1,q2], [CompScienceMeshes.SimplexGraph(1,2)])
translate!(n, point(0,0,20))

X = lagrangec0d1(m, boundary(m))
@test numfunctions(X) == 3

x = refspace(X)
s = chart(m, X.fns[1][1].cellid)
c = neighborhood(s, [1,1]/3)      # get the barycenter of that patch
v = x(c)       # evaluate the Lagrange elements in c, together with their curls

@test (s[3] - s[2]) / (2 * volume(s)) ≈ v[1][2]
@test (s[1] - s[3]) / (2 * volume(s)) ≈ v[2][2]
@test (s[2] - s[1]) / (2 * volume(s)) ≈ v[3][2]

Y = lagrangec0d1(n, boundary(n))
@test numfunctions(Y) == 2

y = refspace(Y)
t = chart(n, Y.fns[1][1].cellid)
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

## tests for the Wilton case
#
m = Mesh([p1,p2,p3], [CompScienceMeshes.SimplexGraph(1,2,3)])
n = boundary(m)#Mesh([p1,p2], [CompScienceMeshes.SimplexGraph(1,2)])

X = lagrangec0d1(m, boundary(m))

Y = lagrangec0d1(n)

# Σ curlfy = 0, assembled matrix nees to have zero row sums
Tyx = assemble(T, Y, X)

@test norm(sum(Tyx,dims=2)) < 1e-12*norm(Tyx)

# finding the threshold distance for our triangle and line combination
# given by xtol2 = 0.2*0.2=0.04 multiplied by 16 in the case κ==0

σ = chart(m, 1)
d_thr = sqrt(0.64*volume(σ))

tol = 1e-6

n_near = translate(n, [0.0,0.0,d_thr*(1-tol)]) # slightly inside threshold for Wilton
n_far = translate(n, [0.0,0.0,d_thr*(1+tol)]) # slightly outside threshold for Wilton

Y_near = lagrangec0d1(n_near)
Y_far = lagrangec0d1(n_far)

Tyx_near = assemble(T, Y_near, X)
Tyx_far = assemble(T, Y_far, X)

# small change in position, but across threshold between two methods, should give similar results
@test norm(Tyx_near - Tyx_far)/norm(Tyx_far) < 1e-5

## the Wilton branch must stay continuous down to d = 0 (test edge ON the triangle)
A0 = Matrix(assemble(T, Y, X)) # Y = lagrangec0d1(boundary(m)), i.e. d = 0
errs = Float64[]
for f in (0.5, 0.1, 0.01, 0.001)
    Yd = lagrangec0d1(translate(n, [0.0,0.0,f*d_thr]))
    push!(errs, norm(Matrix(assemble(T, Yd, X)) - A0)/norm(A0))
end
@test issorted(errs, rev=true)
@test errs[end] < 2e-3

## γ ≠ 0 exercises the regular part (at γ = 0 it vanishes)
nseg = Mesh([p1,p2], [CompScienceMeshes.SimplexGraph(1,2)])
translate!(nseg, point(0,0,0.1*d_thr))
Yseg = lagrangec0d1(nseg, boundary(nseg))

τ = chart(nseg, 1)
tg = normalize(tangents(neighborhood(τ, [0.5]),1))
fv = refspace(X)(neighborhood(σ, [1,1]/3))
# A(γ) - A(0) = -γ·c1 + O(γ²),  c1[i,j] = |σ|·(∫_τ g_i)·⟨t,curl_j⟩/4π,  ∫_τ g_i = |τ|/2
c1 = [ (BEAST.volume(τ)/2)*BEAST.volume(σ)*dot(tg, fv[j].curl)/(4π) for i in 1:2, j in 1:3 ]

A0s = Matrix(assemble(BEAST.NitscheHH3(0.0), Yseg, X))

for (γ, tol) in ((1e-2, 2e-3), (1e-3, 2e-4))
    r = (Matrix(assemble(BEAST.NitscheHH3(γ), Yseg, X)) - A0s) ./ γ
    @test norm(r .+ c1)/norm(c1) < tol
end
