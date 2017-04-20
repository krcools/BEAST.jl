using CompScienceMeshes
using BEAST

κ = ω = 1.0;

Γ = readmesh(Pkg.dir("BEAST", "examples", "sphere.in"))

Xh = raviartthomas(Γ)
Yh = buffachristiansen(Γ)

T = MWSingleLayer3D(κ)
N = NCross()

Txx = assemble(T,Xh,Xh)
@show cond(Txx)
Tyy = assemble(T,Yh,Yh)
@show cond(Tyy)
Nxy = assemble(N,Xh,Yh)
A = Tyy * inv(Nxy) * Txx
@show cond(A)
println("Regression test: cond(A): ",  2.424211950244863)
