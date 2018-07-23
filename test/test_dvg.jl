using BEAST
using CompScienceMeshes

using Test

Γ = meshrectangle(1.0,1.0,1.0)
γ = meshsegment(1.0,1.0,3)
X = raviartthomas(Γ, γ)

x = divergence(X)

@test numfunctions(X) == 2
@test numfunctions(x) == numfunctions(X)
for (_f,_g) in zip(x.fns, X.fns)
    @test length(_f) == length(_g)
    if length(_f) == 1
        c = _f[1].cellid
        _s = chart(Γ, Γ.faces[c])
        @test _f[1].coeff == 1 / volume(_s)
    else
        @test _f[1].coeff + _f[2].coeff == 0
    end
end

edges = skeleton(Γ, 1)
vps = cellpairs(Γ, edges)
rwg = raviartthomas(Γ, vps)
dvg = divergence(rwg)

D = spzeros(Int,numcells(edges), numcells(Γ))
for (_m,_f) in enumerate(dvg.fns)
    for _s in _f
        _e = _s.cellid
        D[_m,_e] = _s.coeff
    end
end

C = connectivity(edges, Γ)
@test size(C) == size(D')


## Test the curl of the linear Lagrange elements
using CompScienceMeshes
using BEAST
using Test

m = meshrectangle(1.0, 1.0, 0.25)
X = duallagrangec0d1(m)
@test numfunctions(X) == numcells(m)

Y = curl(X)
@test numfunctions(X) == numfunctions(Y)
@test isa(Y, BEAST.RTBasis)
@test isa(refspace(Y), BEAST.BEAST.RTRefSpace)

# test the chain property
Z = divergence(Y)
for fn ∈ Z.fns
    total = zero(scalartype(Z))
    for _sh ∈ fn
        total += _sh.coeff
    end
    @test total ≈ 0
end
