using CompScienceMeshes
using BEAST
using Test

for T in [Float32, Float64]
    local m = meshrectangle(T(1.0), T(1.0), T(0.25))

    local X = raviartthomas(m, BEAST.Continuity{:none})
    @test numfunctions(X) == 16*2*3
    @test all(length.(X.fns) .== 1)

    local p = positions(X)

    i = 12
    c = X.fns[i][1].cellid
    local ch = chart(m, c)
    local ctr = cartesian(center(ch))
    @test ctr ≈ p[i]

    for (i,f) in enumerate(X.fns)
        @test f[1].coeff == T(1.0)
    end
end