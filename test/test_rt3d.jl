using CompScienceMeshes
using BEAST

using Test

for T in [Float32, Float64]
    local o, x, y, z = euclidianbasis(3,T)
    tet = simplex(x,y,z,o)

    nbd1 = neighborhood(tet, T.([0,1,1]/3))
    nbd2 = neighborhood(tet, T.([1,0,1]/3))
    nbd3 = neighborhood(tet, T.([1,1,0]/3))
    nbd4 = neighborhood(tet, T.([1,1,1]/3))

    local rs = BEAST.NDLCDRefSpace{T}()
    local fcs = BEAST.faces(tet)
    @test dot(rs(nbd1)[1].value, normal(fcs[1])) > 0
    @test dot(rs(nbd2)[2].value, normal(fcs[2])) > 0
    @test dot(rs(nbd3)[3].value, normal(fcs[3])) > 0
    @test dot(rs(nbd4)[4].value, normal(fcs[4])) > 0

    dot(rs(nbd1)[1].value, normal(fcs[1])) * volume(fcs[1]) ≈ 1
    dot(rs(nbd2)[2].value, normal(fcs[2])) * volume(fcs[2]) ≈ 1
    dot(rs(nbd3)[3].value, normal(fcs[3])) * volume(fcs[2]) ≈ 1
    dot(rs(nbd4)[4].value, normal(fcs[4])) * volume(fcs[4]) ≈ 1

    @test dot(rs(nbd1)[2].value, normal(fcs[1])) * volume(fcs[1]) ≈ 0  atol=eps(T)
    @test dot(rs(nbd1)[3].value, normal(fcs[1])) * volume(fcs[1]) ≈ 0  atol=eps(T)
    @test dot(rs(nbd1)[4].value, normal(fcs[1])) * volume(fcs[1]) ≈ 0  atol=eps(T)

    @test dot(rs(nbd2)[1].value, normal(fcs[2])) * volume(fcs[2]) ≈ 0  atol=eps(T)
    @test dot(rs(nbd2)[3].value, normal(fcs[2])) * volume(fcs[2]) ≈ 0  atol=eps(T)
    @test dot(rs(nbd2)[4].value, normal(fcs[2])) * volume(fcs[2]) ≈ 0  atol=eps(T)

    @test dot(rs(nbd3)[1].value, normal(fcs[3])) * volume(fcs[3]) ≈ 0  atol=eps(T)
    @test dot(rs(nbd3)[2].value, normal(fcs[3])) * volume(fcs[3]) ≈ 0  atol=eps(T)
    @test dot(rs(nbd3)[4].value, normal(fcs[3])) * volume(fcs[3]) ≈ 0  atol=eps(T)

    @test dot(rs(nbd4)[1].value, normal(fcs[4])) * volume(fcs[4]) ≈ 0  atol=eps(T)
    @test dot(rs(nbd4)[2].value, normal(fcs[4])) * volume(fcs[4]) ≈ 0  atol=eps(T)
    @test dot(rs(nbd4)[3].value, normal(fcs[4])) * volume(fcs[4]) ≈ 0  atol=eps(T)


    local ctr = cartesian(center(tet))
    @test dot(cartesian(nbd1)-ctr, normal(fcs[1])) > 0
    @test dot(cartesian(nbd2)-ctr, normal(fcs[2])) > 0
    @test dot(cartesian(nbd3)-ctr, normal(fcs[3])) > 0
    @test dot(cartesian(nbd4)-ctr, normal(fcs[4])) > 0
end