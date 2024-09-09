using Test

using CompScienceMeshes
using BEAST

const e0 = point(0.0,0.0,0.0)
const e1 = point(1.0,0.0,0.0)
const e2 = point(0.0,1.0,0.0)
const e3 = point(0.0,0.0,1.0)

for T in [Float32, Float64]


    p = simplex(
        [
            point(T,0.0,0.0,0.0),
            point(T,1.0,0.0,0.0)
        ],
        Val{1}
    )

    q = simplex(
        [
            point(T,0.0,0.0,0.0),
            point(T,0.5,0.0,0.0)
        ],
        Val{1}
    )

    f = BEAST.LagrangeRefSpace{T,1,2,2}()
    local x = neighborhood(p, T.([0.0]))
    v = f(x)

    @test v[1].value == 0
    @test v[2].value == 1

    @test v[1].derivative == -1
    @test v[2].derivative == +1

    Q = restrict(f, p, q)
    @test Q == [
        T(1.0) T(0.5)
        T(0.0) T(0.5)]


    # Test restriction of RT elements
    # ni, no = 6, 7;
    # ui = transpose([CompScienceMeshes.triangleGaussA[ni] CompScienceMeshes.triangleGaussB[ni] ]);
    # uo = transpose([CompScienceMeshes.triangleGaussA[no] CompScienceMeshes.triangleGaussB[no] ]);
    # wi = triangleGaussW[ni];
    # wo = triangleGaussW[no];
    # universe = Universe(1.0, ui, wi, uo, wo);

    p = simplex([T.(2*e0),T.(2*e1),T.(2*e2)], Val{2})
    x = neighborhood(p,T.([0.5, 0.5]))
    ϕ = BEAST.RTRefSpace{T}()
    v = ϕ(x)

    Q = restrict(ϕ, p, p)
    if T==Float64
    @test Q == Matrix(LinearAlgebra.I, 3, 3)

    q = simplex([T.(e0+e1), T.(2*e1), T.(e1+e2)], Val{2})
    Q = restrict(ϕ, p, q)
    @test Q == [
        2 -1 0
        0  1 0
        0 -1 2] // 4
    end

    # Test restriction of Nedelec elements
    ψ = BEAST.NDRefSpace{T}()
    Q = restrict(ψ, p, p)
    if T==Float64
    @test Q == Matrix(LinearAlgebra.I, 3, 3)

    q = simplex([T.(e0+e1), T.(2*e1), T.(e1+e2)], Val{2})
    Q = restrict(ψ, p, q)
    @test Q == [
        2 -1 0
        0  1 0
        0 -1 2] // 4
    end
end