using CompScienceMeshes, BEAST

Γ = readmesh(Pkg.dir("BEAST", "examples", "sphere.in"))
X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

κ = ω = 1.0; γ = κ*im
T = MWSingleLayer3D(γ)
N = NCross()

Txx = assemble(T,X,X)
Tyy = assemble(T,Y,Y)
Nxy = assemble(N,X,Y)

iNxy = inv(Nxy)
A = iNxy' * Tyy * iNxy * Txx

@show cond(Txx)
@show cond(Tyy)
@show cond(A)

wx, wy, wp = eigvals(Txx), eigvals(Tyy), eigvals(A)
