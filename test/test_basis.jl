## Preamble
using Test

using CompScienceMeshes
using BEAST

## The actual tests
for T in [Float32, Float64]
    κ = ω = T(1.0)

    Γ = meshsegment(T(1.0), T(0.5))
    X = lagrangec0d1(Γ)
    @test numvertices(Γ)-2 == numfunctions(X)

    hypersingular = HyperSingular(κ)
    identityop    = Identity()
    doublelayer   = DoubleLayer(κ)

    @time N = assemble(hypersingular, X, X)
    @time I = assemble(identityop, X, X)

    @test size(N) == (numfunctions(X), numfunctions(X))
    @test size(I) == (numfunctions(X), numfunctions(X))
    @test rank(I) == numfunctions(X)

    @time e = assemble(PlaneWaveNeumann(κ, point(0.0, 1.0)), X)
    @test length(e) == numfunctions(X)

    x1 = N \ e;

    # Testing duallagrangec0d1
    if T == Float64
        Γ1 = meshcircle(T(1.0), T(2.5),3)          # creating a triangle
        Γ2 = barycentric_refinement(Γ1)      # creating the refined mesh

        X1 =duallagrangec0d1(Γ1,Γ2)          # creating the basis functions
        X1 = duallagrangec0d1(Γ1, Γ2, x->false, Val{2})
        @test numcells(Γ1) == numfunctions(X1)  # making sure it is assigned according to the coarse mesh segments
        @test length(X1.fns[1])== 6           # making sure each segment represent 6 shapes inside it
        @test length(X1.fns[numfunctions(X1)])== 6  # making sure the last segment functions contains 6 shapes as well
    end
end

## Test linear Lagrange elements on triangles and the computation of their curl
for T in [Float32, Float64]
    Degr = 1
    Dim1 = 3
    NumF = 3
    f = BEAST.LagrangeRefSpace{T,Degr,Dim1,NumF}()
    sphere = readmesh(joinpath(dirname(@__FILE__),"assets","sphere5.in"),T=T)
    s = chart(sphere, first(cells(sphere)))
    t = neighborhood(s, T.([1,1]/3))
    v = f(t, Val{:withcurl})

    A = volume(s)
    @test v[1][2] == (s[3]-s[2])/2A
    @test v[2][2] == (s[1]-s[3])/2A
    @test v[3][2] == (s[2]-s[1])/2A

    @test v[1][1] ≈ 1/3
    @test v[2][1] ≈ 1/3
    @test v[3][1] ≈ 1/3
end

## Test the construction of continuous linear Lagrange elements on 2D surfaces
using CompScienceMeshes
using BEAST
using Test

for T in [Float32, Float64]
    m = meshrectangle(T(1.0), T(1.0), T(0.5), 3)
    X = lagrangec0d1(m)
    x = refspace(X)

    @test numfunctions(x) == 3

    @test numfunctions(X) == 1
    @test length(X.fns[1]) == 6
end
## test the scalar trace for Lagrange functions
using CompScienceMeshes
using BEAST
using Test
for T in [Float64]
    p1 = point(T,0,0,0)
    p2 = point(T,1,0,0)
    p3 = point(T,0,1,0)
    p4 = point(T,1,1,1)
    c1 = index(1,2,3)

    m = Mesh([p1,p2,p3,p4],[c1])
    b = Mesh([p1,p2], [index(1,2)])
    X = lagrangec0d1(m, boundary(m))
    @test numfunctions(X) == 3

    Y = BEAST.strace(X, b)

    @test numfunctions(X) == 3
    @test numfunctions(Y) == 3

    @test length(Y.fns[1]) == 1
    @test length(Y.fns[2]) == 1
    @test length(Y.fns[3]) == 0

    sh = Y.fns[1][1]; @test (sh.cellid, sh.refid, sh.coeff) == (1, 1, 1.0)
    sh = Y.fns[2][1]; @test (sh.cellid, sh.refid, sh.coeff) == (1, 2, 1.0)

    x = refspace(X)

    cell = chart(m, first(cells(m)))
    face = chart(b, first(cells(b)))
    Q = BEAST.strace(x, cell, 3, face)
    @test Q == [1 0 0; 0 1 0]
end

## test Lagrange construction on Junctions
using CompScienceMeshes
using BEAST
using Test

for T in [Float64]
    m1 = meshrectangle(T(1.0), T(0.5), T(0.5))
    m2 = CompScienceMeshes.rotate(m1, T(0.5π)*[1,0,0])
    m3 = CompScienceMeshes.rotate(m1, T(1.0π)*[1,0,0])
    m = weld(m1, m2, m3)

    X = lagrangec0d1(m)
    x = refspace(X)
    @test numfunctions(X) == 1
    @test length(X.fns[1]) == 9

    p = point(T, 0.5, 0.0, 0.0)
    for _s in X.fns[1]
        _cell = m.faces[_s.cellid]
        patch = chart(m, _cell)
        bary = carttobary(patch, p)
        mp = neighborhood(patch, bary)
    end
end
## Test the dual pieweise constant lagrange elemetns
using CompScienceMeshes
using BEAST
using Test

width, height = 1.0, 1.0
h = 0.5

m = meshrectangle(width, height, h)
b = meshsegment(width, width, 3)

X = duallagrangecxd0(m, b)
@test numfunctions(X) == 2
@test isa(refspace(X), BEAST.LagrangeRefSpace)

@test length(X.fns[1]) == 6
@test length(X.fns[2]) == 12

fine = geometry(X)
for _fn in X.fns
    _n = length(_fn)
    for _sh in _fn
        @test _sh.refid == 1
        cellid = _sh.cellid
        _cell = cells(fine)[cellid]
        ptch = chart(fine, _cell)
        @test _sh.coeff * volume(ptch) ≈ 1/_n
    end
end

## Test the construction of dual piecewise linear, globally continuous elements
using CompScienceMeshes
using BEAST
using Test

m = meshrectangle(1.0, 1.0, 0.25)
j = meshsegment(1.0, 1.0, 3)

X = duallagrangec0d1(m,j)
x = refspace(X)

@test numfunctions(X) == numcells(m)

isonjunction = inclosure_gpredicate(j)
els, ad = BEAST.assemblydata(X)
for _p in 1:numcells(m)
    el = els[_p]
    for r in 1:numfunctions(x)
        vert = el[r]
        isonjunction(vert) || continue
        for (i,w) in ad[_p,r]
            @test w == 0
        end
    end
end



## end of file
