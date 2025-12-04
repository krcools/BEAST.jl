using Test
using LinearAlgebra

using BEAST, CompScienceMeshes, StaticArrays

# Float32 not working since hankelh2 returns always F64
for T in [Float64]

    v1 = SVector(0.0, 0.0)*0.5
    v2 = SVector(1.0, 0.0)*0.5
    v3 = SVector(0.25, 0.0)*0.5
    v4 = SVector(0.50, 0.0)*0.5
    v5 = SVector(0.75, 0.0)*0.5

    shift = SVector(0.5, 0.0)

    v6 = v2 + shift
    v7 = v3 + shift
    v8 = v4 + shift
    v9 = v5 + shift

    verts = [v1, v2, v3, v4, v5, v6, v7, v8, v9]
    faces = [SVector(1, 2, 3, 4, 5), SVector(2, 6, 7, 8, 9)]

    Γto = CompScienceMeshes.CurvilinearMesh(verts, faces, 4)
    Γt = Γto
    Γs = Γt

    Xt = lagrangecxd0(Γt)
    Xs = lagrangecxd0(Γs)

    λ = 10
    k = 2π/λ

    ops = [
        Helmholtz2D.singlelayer(; wavenumber=k)
    ]

    refstrat = BEAST.DoubleNumQStrat(1000, 1001)

    for op in ops
        Sref = assemble(op, Xt, Xs; quadstrat=refstrat)
        ref = Sref[1, 1]

        n = 25

        #Sgl = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumQStrat(n, n+1))
        Sss = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumSauterQstrat(3,3,0,4,n,n))

        #gl = Sgl[1, 1]
        ss = Sss[1, 1]

        #glrel = norm(gl - ref) /  norm(ref)
        ssrel = norm(ss - ref) /  norm(ref)

        @test ssrel < 1e-3
    end

    Xt = lagrangec0d1(Γt)
    Xs = lagrangec0d1(Γs)

    ops = [
        Helmholtz2D.singlelayer(; wavenumber=k)
        Helmholtz2D.hypersingular(; wavenumber=k)
    ]
    𝒮 = 

    refstrat = BEAST.DoubleNumQStrat(1000, 1001)

    for op in ops
        Sref = assemble(op, Xt, Xs; quadstrat=refstrat)
        ref = Sref[1, 1]

        n = 25

        #Sgl = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumQStrat(n, n+1))
        Sss = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumSauterQstrat(3,3,0,4,n,n))

        #gl = Sgl[1, 1]
        ss = Sss[1, 1]

        #glrel = norm(gl - ref) /  norm(ref)
        ssrel = norm(ss - ref) /  norm(ref)

        @test ssrel < 1e-3
    end

    ## Test accuracy of vertex integrations

    v1 = SVector(0.0, 0.0)
    v2 = SVector(1.0, 0.0)
    v3 = SVector(0.25, 0.0)
    v4 = SVector(0.50, 0.0)
    v5 = SVector(0.75, 0.0)
    v6 = SVector(0.0, 1.0)
    v7 = SVector(0.0, 0.25)
    v8 = SVector(0.0, 0.50)
    v9 = SVector(0.0, 0.75)

    verts = [v1, v2, v3, v4, v5, v6, v7, v8, v9]
    facest = [SVector(1, 2, 3, 4, 5)]
    facess = [SVector(1, 6, 7, 8, 9)]

    Γt = CompScienceMeshes.CurvilinearMesh(verts, facest, 4)
    Γs = CompScienceMeshes.CurvilinearMesh(verts, facess, 4)

    Xt = lagrangecx(Γt, order=3)
    Xs = lagrangecx(Γs, order=3)

    λ = 10
    k = 2π/λ

    ops = [
        Helmholtz2D.singlelayer(; wavenumber=k)
        Helmholtz2D.doublelayer(; wavenumber=k)
        Helmholtz2D.doublelayer_transposed(; wavenumber=k)
    ]

    refstrat = BEAST.DoubleNumQStrat(1000, 1001)

    for op in ops
        Sref = assemble(op, Xt, Xs; quadstrat=refstrat)
        ref = Sref[1, 1]

        n = 25

        #Sgl = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumQStrat(n, n+1))
        Sss = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumSauterQstrat(3,3,0,4,n,n))

        #gl = Sgl[1, 1]
        ss = Sss[1, 1]

        #glrel = norm(gl - ref) /  norm(ref)
        ssrel = norm(ss - ref) /  norm(ref)

        @test ssrel < 1e-4

        Ssrel = norm(Sss - Sref)/norm(Sref)

        @test Ssrel < 1e-4
    end


##
    facess = [SVector(6, 1, 9, 8, 7)]

    Γs = CompScienceMeshes.CurvilinearMesh(verts, facess, 4)

    Xt = lagrangecx(Γt, order=3)
    Xs = lagrangecx(Γs, order=3)

    λ = 10
    k = 2π/λ

    ops = [
        Helmholtz2D.singlelayer(; wavenumber=k)
        Helmholtz2D.doublelayer(; wavenumber=k)
        Helmholtz2D.doublelayer_transposed(; wavenumber=k)
    ]

    refstrat = BEAST.DoubleNumQStrat(1000, 1001)

    for op in ops
        Sref = assemble(op, Xt, Xs; quadstrat=refstrat)
        ref = Sref[1, 1]

        n = 25

        #Sgl = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumQStrat(n, n+1))
        Sss = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumSauterQstrat(3,3,0,4,n,n))

        #gl = Sgl[1, 1]
        ss = Sss[1, 1]

        #glrel = norm(gl - ref) /  norm(ref)
        ssrel = norm(ss - ref) /  norm(ref)

        @test ssrel < 1e-4

        Ssrel = norm(Sss - Sref)/norm(Sref)
        @test Ssrel < 1e-4
    end
end