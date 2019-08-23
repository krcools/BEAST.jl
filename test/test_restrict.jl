using Test

using CompScienceMeshes
using BEAST

const e0 = point(0.0,0.0,0.0)
const e1 = point(1.0,0.0,0.0)
const e2 = point(0.0,1.0,0.0)
const e3 = point(0.0,0.0,1.0)

p = simplex(
    [
        point(0.0,0.0,0.0),
        point(1.0,0.0,0.0)
    ],
    Val{1}
)

q = simplex(
    [
        point(0.0,0.0,0.0),
        point(0.5,0.0,0.0)
    ],
    Val{1}
)

f = BEAST.LagrangeRefSpace{Float64,1,2,2}()
x = neighborhood(p, [0.0])
v = f(x)

@test v[1].value == 0
@test v[2].value == 1

@test v[1].derivative == -1
@test v[2].derivative == +1

Q = restrict(f, p, q)
@test Q == [
    1.0 0.5
    0.0 0.5]


# Test restriction of RT elements
ni, no = 6, 7;
ui = transpose([triangleGaussA[ni] triangleGaussB[ni] ]);
uo = transpose([triangleGaussA[no] triangleGaussB[no] ]);
wi = triangleGaussW[ni];
wo = triangleGaussW[no];
# universe = Universe(1.0, ui, wi, uo, wo);

p = simplex([2*e0,2*e1,2*e2], Val{2})
x = neighborhood(p,[0.5, 0.5])
ϕ = BEAST.RTRefSpace{Float64}()
v = ϕ(x)

Q = restrict(ϕ, p, p)
@test Q == Matrix(LinearAlgebra.I, 3, 3)

q = simplex([e0+e1, 2*e1, e1+e2], Val{2})
Q = restrict(ϕ, p, q)
@test Q == [
    2 -1 0
    0  1 0
    0 -1 2] // 4
