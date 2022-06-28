using CompScienceMeshes
using BEAST
using Test
using LinearAlgebra

for T in [Float32, Float64]
    local o, x, y, z = CompScienceMeshes.euclidianbasis(3,T)
    tet = simplex(x,y,z,o)

    rs = BEAST.NDLCDRefSpace{T}()
    Q = BEAST.restrict(rs, tet, tet)
    @test Q ≈ Matrix(T(1.0)LinearAlgebra.I, 4, 4)

    rs = BEAST.NDLCCRefSpace{T}()
    Q = BEAST.restrict(rs, tet, tet)
    @test Q ≈ Matrix(T(1.0)LinearAlgebra.I, 6, 6)

    c = cartesian(center(tet))
    smalltet = simplex(x,y,z,c)
    p_smalltet = center(smalltet)
    p_tet = neighborhood(tet, carttobary(tet, cartesian(p_smalltet)))
    @show cartesian(p_smalltet)
    @show cartesian(p_tet)
    @assert cartesian(p_tet) ≈ cartesian(p_smalltet)
    @assert volume(tet) / volume(smalltet) ≈ 4

    Q = BEAST.restrict(rs, tet, smalltet)
    Fp = rs(p_tet)
    fp = rs(p_smalltet)

    for j in axes(Q,1)
        x = Fp[j].value
        y = sum(Q[j,i]*fp[i].value for i in axes(Q,2))
        @show x
        @show y
        @test x ≈ y
    end
end