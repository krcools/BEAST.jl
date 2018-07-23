using Test
#using LinearForms
using CompScienceMeshes
using BEAST

x = point(1.0, 0.0, 0.0)
y = point(0.0, 1.0, 0.0)
z = point(0.0, 0.0, 1.0)

κ = 1.0
S = SingleLayerTrace(im*κ)

h = 0.25
Γ = meshrectangle(1.0,1.0,h)
γ = meshsegment(1.0,1.0,3)

X = raviartthomas(Γ, γ)

x = divergence(X)
y = ntrace(X,γ)
Z = assemble(S,y,x)

# test for the correct sparsity pattern
# I, J, V = findall(!iszero, Z)
Q = findall(!iszero, Z)
I = getindex.(Q,1)
J = getindex.(Q,2)
@test length(unique(I)) == 4
@test length(unique(J)) == 44


## test the acutal value of the penalty terms
using CompScienceMeshes
using BEAST
using Test

p1 = point(0,0,0)
p2 = point(1,0,0)
p3 = point(0,1,0)

q1 = point(0,0,0)
q2 = point(1,0,0)

m = Mesh([p1,p2,p3],[index(1,2,3)])
n = Mesh([q1,q2], [index(1,2)])
translate!(n, point(0,0,20))

X = lagrangecxd0(m)
Y = lagrangecxd0(n)

x = refspace(X)
y = refspace(Y)

N = BEAST.SingleLayerTrace(0.0)
Nyx = assemble(N,Y,X)

@test size(Nyx) == (1,1)

sx = chart(m, first(cells(m)))
sy = chart(n, first(cells(n)))

cx = neighborhood(sx, [1,1]/3)
cy = neighborhood(sy, [1]/2)

vx = x(cx)
vy = y(cy)

R = norm(cartesian(cx) - cartesian(cy))
estimate = volume(sx) * volume(sy) * vx[1][1] * vy[1][1] / (4π*R)
actual = Nyx[1,1]
@test norm(estimate - actual) / norm(actual) < 5e-4

# test that the trace and divergence work as advetised
X = raviartthomas(m,boundary(m))

D = divergence(X)
Y = ntrace(X,boundary(m))

@test numfunctions(D) == 3
@test isa(refspace(X), BEAST.RTRefSpace)
for _f in D.fns
    @test length(_f) == 1
    @test _f[1].cellid == 1
    @test _f[1].refid == 1
    @test _f[1].coeff ≈ (1 / volume(sx))
end

Σ = geometry(Y)
@test numfunctions(Y) == 3
@test isa(refspace(Y), BEAST.LagrangeRefSpace)
for _f in Y.fns
    @test length(_f) == 1
    @test 0 < _f[1].cellid < 4
    seg = chart(Σ, Σ.faces[_f[1].cellid])
    @test _f[1].refid == 1
    @test _f[1].coeff ≈ (1 / volume(seg))
end
