using Test

using BEAST
using CompScienceMeshes

chart1 = simplex(
    point(1,0,0),
    point(0,1,0),
    point(0,0,0))

chart2 = simplex(
    point(1/2,0,0),
    point(0,1/2,0),
    point(1,1,0))

X = BEAST.RTRefSpace{Float64}()
@time Q1 = BEAST.restrict(X, chart1, chart2)
@time Q2 = BEAST.interpolate(X, chart2, X, chart1)
@test Q1 ≈ Q2

constant_vector_field = point(1,2,0)
Q3 = BEAST.interpolate(X, chart2) do p
    return [constant_vector_field]
end

ctr = center(chart2)
vals = [f.value for f in X(ctr)]
itpol = sum(w*val for (w,val) in zip(Q3,vals))
@test itpol ≈ constant_vector_field

# using TestItems
@testitem "restrict RT0" begin
    using CompScienceMeshes

    ref_vertices = [
        point(1,0),
        point(0,1),
        point(0,0),
    ]
    vertices = [
        point(1,0,0),
        point(0,1,0),
        point(0,0,0),
    ]
    chart1 = simplex(vertices...)
    for I in BEAST._dof_perms_rt
        chart2 = simplex(
            chart1.vertices[I[1]],
            chart1.vertices[I[2]],
            chart1.vertices[I[3]],)
        chart2tochart1 = CompScienceMeshes.simplex(ref_vertices[collect(I)]...)
        rs = BEAST.RTRefSpace{Float64}()
        Q1 = BEAST.dof_perm_matrix(rs, I)
        Q2 = BEAST.restrict(rs, chart1, chart2)
        Q3 = BEAST.restrict(rs, chart1, chart2, chart2tochart1)
        @test Q1 ≈ Q2
        @test Q1 ≈ Q3
    end
end

@testitem "restrict Lagrange on a refined cell" begin
    using CompScienceMeshes
    using LinearAlgebra

    coarse = simplex(point(0,0,0), point(1,0,0), point(0,1,0))
    # one child of the barycentric refinement of coarse:
    child = simplex(point(1,0,0),point(0.5,0.5,0), point(1/3,1/3,0))

    tobary(chart, ch) = simplex(map(v -> carttobary(chart,v), ch.vertices))

    # the three routes to a restriction matrix must agree on a genuine
    # restriciton; self_restriction alone passes under any node ordering
    for (degree, n) in ((1,3),(2,6))
        rs = BEAST.LagrangeRefSpace{Float64, degree, 3, n}()
        Q1 = BEAST.restrict(rs, coarse, child)
        Q2 = BEAST.restrict(rs, coarse, child, tobary(coarse, child))
        Q3 = zeros(n, n)
        BEAST.restrict!(Q3, rs, coarse, child, tobary(coarse, child))
        @test Q1 ≈ Q2
        @test Q1 ≈ Q3
    end

    for (degree, n) in ((0,1),(1,3),(2,6),(3,10))
        rs = BEAST.LagrangeRefSpace{Float64, degree, 3, n}()
        Q = zeros(n, n)
        BEAST.restrict!(Q, rs, coarse, coarse, tobary(coarse, coarse))
        @test Q ≈ I
    end
end
