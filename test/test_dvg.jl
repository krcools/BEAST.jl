using BEAST
using CompScienceMeshes

using Base.Test

Γ = meshrectangle(1.0,1.0,1.0)
γ = meshsegment(1.0,1.0,3)
X = raviartthomas(Γ, γ)

x = divergence(X)

@test numfunctions(X) == 2
@test numfunctions(x) == numfunctions(X)
for (f,g) in zip(x.fns, X.fns)
    @test length(f) == length(g)
    if length(f) == 1
        c = f[1].cellid
        #s = simplex(cellvertices(Γ,c))
        s = simplex(vertices(Γ, Γ.faces[c]))
        @test f[1].coeff == 1 / volume(s)
    else
        @test f[1].coeff + f[2].coeff == 0
    end
end

edges = skeleton(Γ, 1)
vps = cellpairs(Γ, edges)
rwg = raviartthomas(Γ, vps)
dvg = divergence(rwg)

D = spzeros(Int,numcells(edges), numcells(Γ))
for (m,f) in enumerate(dvg.fns)
    for s in f
        e = s.cellid
        D[m,e] = s.coeff
    end
end

C = connectivity(edges, Γ)
@test size(C) == size(D')


## Test the curl of the linear Lagrange elements
using CompScienceMeshes
using BEAST
using Base.Test

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
    for sh ∈ fn
        total += sh.coeff
    end
    @test total ≈ 0
end
