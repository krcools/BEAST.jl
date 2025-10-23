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