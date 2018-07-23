using CompScienceMeshes
using BEAST

using Test
using StaticArrays

T = Float64
P = SVector{3,T}

tol = eps(T) * 10^3

function shapevals(ϕ, ts)

    numpoints = length(ts)
    numshapes = numfunctions(ϕ)

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

const third = one(T)/3
m = meshrectangle(1.0,1.0,0.5);
rwg = raviartthomas(m)
ref = refspace(rwg)

# vrts = vertices(m, first(cells(m)))
# ptch = simplex(vrts[1], vrts[2], vrts[3])
ptch = chart(m, first(cells(m)))

ctrd = neighborhood(ptch, [1,1]/3)
vals = shapevals(ref, [ctrd])

# test edge detection
edges = skeleton(m, 1)
vp = cellpairs(m,edges)

# test function evaluation
v = [
    point(0.0, 0.0, 0.0),
    point(2.0, 0.0, 0.0),
    point(0.0, 3.0, 0.0)]

cell = simplex(v[1], v[2], v[3])


mp = neighborhood(cell,[third, third]);
ϕ = BEAST.RTRefSpace{Float64}()
fx = ϕ(mp)[2][1]

A  = volume(cell)
r = cartesian(mp)
gx = (r - v[2]) / 2A
@test norm(fx - gx) < tol

@test norm((r - v[1])/2A -  ϕ(mp)[1][1]) < tol
@test norm((r - v[2])/2A -  ϕ(mp)[2][1]) < tol
@test norm((r - v[3])/2A -  ϕ(mp)[3][1]) < tol

# repeat the test but now on a random cell
randpoint() = point(2*rand()-1, 2*rand()-1, 2*rand()-1)
v = [randpoint(), randpoint(), randpoint()]
cell = simplex(v[1], v[2], v[3])
mp = neighborhood(cell,[third, third]);

fx = ϕ(mp)[2][1]
#fx = evalfun(rs, mp)

A = volume(cell)
r = cartesian(mp)
gx = (r-v[2]) / 2A
@test norm(fx - gx) < tol

@test norm((r - v[1])/2A -  ϕ(mp)[1][1]) < tol
@test norm((r - v[2])/2A -  ϕ(mp)[2][1]) < tol
@test norm((r - v[3])/2A -  ϕ(mp)[3][1]) < tol
