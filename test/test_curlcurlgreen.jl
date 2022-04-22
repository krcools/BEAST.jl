using Test
using LinearAlgebra
using BEAST
using CompScienceMeshes
using StaticArrays
for T in [Float32, Float64]
    local y = point(T,2,0,0)
    local j = ẑ
    local κ = T(1.0)

    gj(x) = exp(-im*κ*norm(x-y))/(4*π*norm(x-y)) * j
    ccg = BEAST.CurlCurlGreen(κ, j, y)

    function curlh(f,x,h)
        e1 = point(T,1,0,0)
        e2 = point(T,0,1,0)
        e3 = point(T,0,0,1)

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

    local h = T(0.01)
    local x = point(T,1,1,1)
    a = ccg(x)
    local b = ccgh(x)
    @show norm(x-y)
    T == Float64 ? atol=1e-5 : atol=1e-4
    @test norm(a-b) < atol

    cccg = curl(ccg)
    cccgh(x) = curlh(ccg,x,h)

    # @show a = cccg(x)
    # @show b = cccgh(x)
    @test norm(a-b) < atol
end