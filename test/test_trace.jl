using CompScienceMeshes
using BEAST
using Test

m = meshrectangle(1.0, 1.0, 0.5)
b = meshsegment(1.0, 0.3, 3)

edges = skeleton(m, 1)

vps = cellpairs(m, edges)
rwg = raviartthomas(m, vps)

nt = ntrace(rwg, b)

#nt = BEAST.ntrace2(rwg, b)

n = 0
Σ = geometry(nt)
on_bnd = overlap_gpredicate(b)
for i in eachindex(rwg.fns)
    length(nt.fns[i]) == 1 || continue
    global n += 1
    c = nt.fns[i][1].cellid
    edge = chart(Σ, Σ.faces[c])
    @test on_bnd(edge)
    @test nt.fns[i][1].refid == 1
    @test nt.fns[i][1].coeff == 2.0
end

@test n == 2

## test the scalar trace of a lgrange basis
using CompScienceMeshes
using BEAST
using Test
m = meshrectangle(1.0, 1.0, 0.5)
b = meshsegment(1.0, 0.3, 3)

X = lagrangec0d1(m,b)

@test numfunctions(X) == 2
@test length(X.fns[1]) == 3
@test length(X.fns[2]) == 6

Y = strace(X,b)
#A = findall(length(f) for f in Y.fns)
A = findall(length.(Y.fns) .!= 0)

@test length(A) == 1
@test length(Y.fns[1]) == 2
@test Y.fns[1][1].coeff ≈ 1
@test Y.fns[1][2].coeff ≈ 1
