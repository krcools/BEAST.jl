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

using TestItems
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


