using CompScienceMeshes
using BEAST
using LinearAlgebra
using Test
using StaticArrays

for T in [Float32, Float64]
    local o, x, y, z = euclidianbasis(3, T)

    p1 = 2x
    p2 = y
    p3 = 3z
    p4 = o
    tet = simplex(p1, p2, p3, p4)

    q = abs(CompScienceMeshes.relorientation([1,2,3],[1,2,3,4]))
    tri = simplex(p1, p2, p3)
    local n = normal(tri)

    i = 3
    local j = 3
    local edg = simplex(p1, p2)


    x = BEAST.NDLCDRefSpace{T}()
    y = BEAST.NDRefSpace{T}()

    local p = neighborhood(tet, T.([0.5, 0.5, 0.0]))
    r = neighborhood(tri, T.([0.5, 0.5]))

    @test carttobary(edg, cartesian(p)) ≈ T.([0.5])
    @test carttobary(edg, cartesian(r)) ≈ T.([0.5])

    xp = x(p)[j].value
    yr = y(r)[i].value

    local a, b = extrema((n × xp) ./ yr)
    @test a ≈ b

    tgt = p2 - p1
    a2 = dot(n × xp, tgt) / dot(yr, tgt)
    @test a ≈ a2

    z = a * yr
    @test z ≈ n × xp

    volume(simplex(p1,p2,p4))
    volume(simplex(p1,p2))

    Q = BEAST.ttrace(x, tet, q, tri)

    p = neighborhood(tet, T.([1/3, 1/3, 1/3]))
    r = neighborhood(tri, T.([1/3, 1/3]))
    @test cartesian(p) ≈ cartesian(r)

    xp = n × x(p)[j].value
    yr = y(r)[i].value
    @test xp ≈ yr * Q[i,j]


    x = BEAST.NDLCCRefSpace{T}()
    y = BEAST.RTRefSpace{T}()
    for q in 1:4
        q = 3
        fc = BEAST.faces(tet)[q]
        Q = BEAST.ttrace(x, tet, q, fc)

        nbdi = CompScienceMeshes.center(fc)
        nbdj = neighborhood(tet, carttobary(tet, cartesian(nbdi)))

        xvals = x(nbdj)
        yvals = y(nbdi)

        for j in 1:6
            trc = sum(Q[i,j]*yvals[i].value for i in 1:3)
            @test isapprox(trc, normal(fc) × xvals[j].value, atol=1e-4)
        end
    end

    # test the case where intrinsic and extrinsic orientations differ
    o, x, y, z = euclidianbasis(3, T)
    fc = simplex(z,o,y)
    x = BEAST.NDLCCRefSpace{T}()
    y = BEAST.RTRefSpace{T}()
    Q = BEAST.ttrace(x, tet, 3000, fc)

    nbdi = CompScienceMeshes.center(fc)
    nbdj = neighborhood(tet, carttobary(tet, cartesian(nbdi)))

    @test cartesian(nbdi) ≈ cartesian(nbdj)

    xvals = x(nbdj)
    yvals = y(nbdi)

    for j in 1:6
        trc = sum(Q[i,j]*yvals[i].value for i in 1:3)
        @test isapprox(trc, -normal(fc) × xvals[j].value, atol=1e-4)
    end


    o, x, y, z = euclidianbasis(3)
    local m = Mesh([x,y,z,o], [CompScienceMeshes.SimplexGraph(1,2,3,4)])

    m1 = skeleton(m,1)
    local X = BEAST.nedelecc3d(m, m1)
    @test numfunctions(X) == 6

    # m2 = skeleton(m,2)
    m2 = boundary(m)
    local Y = BEAST.ttrace(X,m2)
    @test numfunctions(Y) == 6

    pa = Y.fns[1][1].cellid
    pb = Y.fns[1][2].cellid

    tria = chart(m2, pa)
    trib = chart(m2, pb)

    ra = Y.fns[1][1].refid
    rb = Y.fns[1][2].refid

    ctra = cartesian(CompScienceMeshes.center(CompScienceMeshes.edges(tria)[ra]))
    ctrb = cartesian(CompScienceMeshes.center(CompScienceMeshes.edges(trib)[rb]))
    @test ctra ≈ ctrb
    # tri1 = chart(m, cells(m)[Y.fns[1][1].cellid])

    # ctr1 = cartesian(center(CompScienceMeshes.edges(tri)[]))
end