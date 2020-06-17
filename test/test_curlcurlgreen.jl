using Test
using LinearAlgebra
using BEAST
using CompScienceMeshes
using StaticArrays

y = point(2,0,0)
j = ẑ
κ = 1.0

gj(x) = exp(-im*κ*norm(x-y))/(4*π*norm(x-y)) * j
ccg = BEAST.CurlCurlGreen(κ, j, y)

function curlh(f,x,h)
    e1 = point(1,0,0)
    e2 = point(0,1,0)
    e3 = point(0,0,1)

    d1f = (f(x+h*e1) - f(x-h*e1))/(2*h)
    d2f = (f(x+h*e2) - f(x-h*e2))/(2*h)
    d3f = (f(x+h*e3) - f(x-h*e3))/(2*h)

    return @SVector[
        d2f[3] - d3f[2],
        d3f[1] - d1f[3],
        d1f[2] - d2f[1]]
end

cgh(x) = curlh(gj,x,h)
ccgh(x) = curlh(cgh,x,h)

h = 0.01
x = point(1,1,1)
a = ccg(x)
b = ccgh(x)
@show norm(x-y)
@test norm(a-b) < 1e-5

cccg = curl(ccg)
cccgh(x) = curlh(ccg,x,h)

# @show a = cccg(x)
# @show b = cccgh(x)
@test norm(a-b) < 1e-5
