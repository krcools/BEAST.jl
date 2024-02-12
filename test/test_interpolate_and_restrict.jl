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