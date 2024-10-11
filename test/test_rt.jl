using CompScienceMeshes
using BEAST

using Test
using StaticArrays

const third = one(Float64)/3
for T in [Float32, Float64]
P = SVector{3,T}

tol = eps(T) * 10^3

function shapevals(ϕ, ts)

    numpoints = length(ts)
    dom = CompScienceMeshes.domain(CompScienceMeshes.chart(ts[1]))
    numshapes = numfunctions(ϕ, dom)

    y = Array{P}(undef, numshapes, numpoints)
    for i in 1 : numpoints
        t = ts[i]
        u = ϕ(t)
        for j in 1 : numshapes
            y[j,i] = u[j][1]
        end
    end

    return y
end


local m = meshrectangle(T(1.0),T(1.0),T(0.5));
rwg = raviartthomas(m)
ref = refspace(rwg)

# vrts = vertices(m, first(cells(m)))
# ptch = simplex(vrts[1], vrts[2], vrts[3])
ptch = chart(m, first(m))

ctrd = neighborhood(ptch, T.([1,1]/3))
local vals = shapevals(ref, [ctrd])

# test edge detection
edges = skeleton(m, 1)
vp = cellpairs(m,edges)

# test function evaluation
v = [
    point(T, 0.0, 0.0, 0.0),
    point(T, 2.0, 0.0, 0.0),
    point(T, 0.0, 3.0, 0.0)]

cell = simplex(v[1], v[2], v[3])


mp = neighborhood(cell,[T(third), T(third)]);
ϕ = BEAST.RTRefSpace{T}()
fx = ϕ(mp)[2][1]

A  = volume(cell)
r = cartesian(mp)
gx = (r - v[2]) / 2A
@test norm(fx - gx) < tol

@test norm((r - v[1])/2A -  ϕ(mp)[1][1]) < tol
@test norm((r - v[2])/2A -  ϕ(mp)[2][1]) < tol
@test norm((r - v[3])/2A -  ϕ(mp)[3][1]) < tol

# repeat the test but now on a random cell
randpoint() = point(2*rand(T)-1, 2*rand(T)-1, 2*rand(T)-1)
v = [randpoint(), randpoint(), randpoint()]
cell = simplex(v[1], v[2], v[3])
mp = neighborhood(cell,[T(third), T(third)]);

fx = ϕ(mp)[2][1]
#fx = evalfun(rs, mp)

A = volume(cell)
r = cartesian(mp)
gx = (r-v[2]) / 2A
@test norm(fx - gx) < tol

@test norm((r - v[1])/2A -  ϕ(mp)[1][1]) < tol
@test norm((r - v[2])/2A -  ϕ(mp)[2][1]) < tol
@test norm((r - v[3])/2A -  ϕ(mp)[3][1]) < tol
end