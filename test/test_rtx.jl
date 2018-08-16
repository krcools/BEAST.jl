using CompScienceMeshes
using BEAST
using Test

m = meshrectangle(1.0, 1.0, 0.25)

X = raviartthomas(m, BEAST.Continuity{:none})
@test numfunctions(X) == 16*2*3
@test all(length.(X.fns) .== 1)

p = positions(X)

i = 12
c = X.fns[i][1].cellid
ch = chart(m, cells(m)[c])
ctr = cartesian(center(ch))
@test ctr ≈ p[i]

for (i,f) in enumerate(X.fns)
    @test f[1].coeff == 1.0
end
