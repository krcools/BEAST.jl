using CompScienceMeshes, BEAST

Γ = readmesh(Pkg.dir("BEAST", "examples", "sphere.in"))
X, Y = raviartthomas(Γ), buffachristiansen(Γ)

κ = ω = 1.0; γ = κ*im
T = MWSingleLayer3D(γ)
N = NCross()

Txx = sdata(assemble(T,X,X))
Tyy = sdata(assemble(T,Y,Y))
Nxy = sdata(assemble(N,X,Y))

iNxy = inv(Nxy)
A = iNxy' * Tyy * iNxy * Txx

@show cond(Txx)
@show cond(Tyy)
@show cond(A)

wx, wy, wp = eigvals(Txx), eigvals(Tyy), eigvals(A)
