using Test

using CompScienceMeshes
using BEAST

const e0 = Point(0.0,0.0,0.0)
const e1 = Point(1.0,0.0,0.0)
const e2 = Point(0.0,1.0,0.0)
const e3 = Point(0.0,0.0,1.0)

p = patch(
    [
        Point(0.0,0.0,0.0),
        Point(1.0,0.0,0.0)
    ],
    Val{1}
)

q = patch(
    [
        Point(0.0,0.0,0.0),
        Point(0.5,0.0,0.0)
    ],
    Val{1}
)

f = LagrangeRefSpace{Float64,1}()
x = neighborhood(p, [0.0])
v = f(x)
@test v == Vec(0.0, 1.0)

Q = restrict(f, p, q)
@test Q == [
    1.0 0.5
    0.0 0.5
]


# Test restriction of RT elements
ni, no = 6, 7;
ui = transpose([triangleGaussA[ni] triangleGaussB[ni] ]);
uo = transpose([triangleGaussA[no] triangleGaussB[no] ]);
wi = triangleGaussW[ni];
wo = triangleGaussW[no];
universe = Universe(1.0, ui, wi, uo, wo);

p = patch([2*e0,2*e1,2*e2], Val{2})
x = neighborhood(p,[0.5, 0.5])
ϕ = RTRefSpace{Float64}()
v = ϕ(x)

Q = restrict(ϕ, p, p)
@test Q == eye(3)

q = patch([e0+e1, 2*e1, e1+e2], Val{2})
Q = restrict(ϕ, p, q)
@test Q == [
    2 -1 0
    0  1 0
    0 -1 2] // 4
