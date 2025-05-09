@testitem "excitation: TD Helmholtz 3D" begin
    # using BEAST
    using CompScienceMeshes
    # using Test
    for T in [Float32, Float64]
        dir = point(T,0,0,1)
        width = T(1.0)
        delay = T(1.5)
        scaling = T(1.0)
        sig = creategaussian(width, delay, scaling)

        pw = BEAST.planewave(dir, T(1.0), sig)

        CompScienceMeshes.cartesian(p::typeof(dir)) = p
        CompScienceMeshes.cartesian(p::Number) = p
        x = point(T,1,0,0)
        t = T(1.0)
        @test pw(x,t) ≈ sig(t-dot(dir,x))

        pw2 = BEAST.gradient(pw)
        @test pw2.direction ≈ dir
        @test pw2.polarisation ≈ -dir
        @test pw2.speedoflight ≈ pw.speed_of_light

        dsig = derive(sig)
        @test pw2(x,t) ≈ -dir*dsig(t-dot(dir,x))

        trc = dot(BEAST.n,pw2)
        ch = simplex(
            point(T,0,0,0),
            point(T,1,0,0),
            point(T,0,1,0))
        ctr = center(ch)
        val = trc(ctr,t)
        @test val ≈ -dsig(t-dot(dir,x))
    end
end