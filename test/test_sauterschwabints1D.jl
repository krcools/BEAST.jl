using Test
using LinearAlgebra

using BEAST, CompScienceMeshes, StaticArrays


# Float32 not working since hankelh2 returns always F64
for T in [Float64]
    T = Float64
    Γto = meshsegment(T(1.0), T(0.5))
    Γt = Γto
    Γs = Γt

    Xt = lagrangecxd0(Γt)
    Xs = lagrangecxd0(Γs)

    λ = 10
    k = 2π/λ

    ops = [
        Helmholtz2D.singlelayer(; wavenumber=k)
    ]

    refstrat = BEAST.DoubleNumSauterQstrat(3,3,0,4,30,30)

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

        @test ssrel < 1e-14
    end

    Xt = lagrangec0d1(Γt)
    Xs = lagrangec0d1(Γs)

    ops = [
        Helmholtz2D.singlelayer(; wavenumber=k)
        Helmholtz2D.hypersingular(; wavenumber=k)
    ]
    𝒮 = 

    refstrat = BEAST.DoubleNumSauterQstrat(3,3,0,4,30,30)

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

        @test ssrel < 1e-14
    end

    ## Test accuracy of vertex integrations

    T = Float64
    Γto = meshsegment(T(1.0), T(1.0))
    Γt = Γto
    Γs = translate(Γto, SVector(1.0, 0.0))

    Γs = Mesh([SVector(0.0, 0.0), SVector(0.0, 1.0)], [SVector(1, 2)])
    Xt = lagrangecxd0(Γt)
    Xs = lagrangecxd0(Γs)

    λ = 10
    k = 2π/λ

    ops = [
        Helmholtz2D.singlelayer(; wavenumber=k)
        Helmholtz2D.doublelayer(; wavenumber=k)
        Helmholtz2D.doublelayer_transposed(; wavenumber=k)
    ]
    𝒮 = 

    refstrat = BEAST.DoubleNumSauterQstrat(3,3,0,4,30,30)

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

        @test ssrel < 1e-14
    end
end