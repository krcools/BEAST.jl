using Test

using CompScienceMeshes
using BEAST

using LinearAlgebra


Γ = readmesh(joinpath(@__DIR__, "assets", "sphere2.in"))
X = raviartthomas(Γ);

κ, η = 1.0, 1.0;
t = Maxwell3D.singlelayer(wavenumber=κ);
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ);
e = (n × E) × n;


A = assemble(t,X,X)
b = assemble(e,X)

u_direct = A \ b

res_direct = norm(A*u_direct-b)/norm(b)
@test res_direct < eps()*100

for rtol in [1e-6,1e-8,1e-10]

    u = BEAST.GMRESSolver(A; reltol=rtol) * b


    res = norm(A*u-b)/norm(b)

    diff = norm(u-u_direct)/norm(u)

    @test res < rtol 
    @test diff < rtol*100
end


@testitem "GMRESSolver: left/right_preconditioner kwarg" begin
    using LinearAlgebra
     A = [
        0.79569   0.484796  0.68263   0.741895  0.936866
        0.222479  0.215901  0.441589  0.55834   0.572374
        0.165348  0.251567  0.789082  0.456141  0.225406
        0.442217  0.25817   0.328375  0.922877  0.476834
        0.154665  0.980761  0.58988   0.881435  0.014069
     ]

     b = [
        0.7261143955334382
        0.14878384999306027
        0.33554623017710794
        0.15630909408720672
        0.06524871057556192
     ]

    Ai1 = inv(A)
    Ai2 = BEAST.GMRESSolver(A; reltol=1e-8)

     u0 = A \ b
     u1, ch1 = solve(Ai2, b)
     u2, ch2 = solve(Ai2, b; left_preconditioner=Ai1)

    @test norm(u1 - u0) / norm(u0) < 1e-7
    @test norm(u2 - u0) / norm(u0) < 1e-7

    @test ch2.iters == 1

end