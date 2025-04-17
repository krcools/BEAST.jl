using CompScienceMeshes, BEAST
using LinearAlgebra

Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
X = raviartthomas(Γ)
Y = BEAST.buffachristiansen(Γ)

κ = 1.0;
T = Maxwell3D.singlelayer(wavenumber=κ)
N = NCross()

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

bx = assemble(e, X)

A = assemble(T,X,X); println("primal discretisation assembled.")

Tyy = assemble(T,Y,Y); println("dual discretisation assembled.")
Nxy = Matrix(assemble(N,X,Y)); println("duality form assembled.")
iNxy = inv(Nxy); println("duality form inverted.")
P = iNxy' * Tyy * iNxy

iT1 = BEAST.GMRESSolver(A; restart=1_500, abstol=1e-8, reltol=1e-8, maxiter=1_500)
iT2 = BEAST.GMRESSolver(A; restart=1_500, abstol=1e-8, reltol=1e-8, maxiter=1_500, left_preconditioner=P)

u1, ch1 = BEAST.solve(iT1, bx)
u2, ch2 = BEAST.solve(iT2, bx)

using Plots
Plots.plot(title="log10 residual error vs iteration count")
Plots.plot!(log10.(ch1[:resnorm]), label="A=b", l=2)
Plots.plot!(log10.(ch2[:resnorm]), label="P*A=P*b", l=2)

ss1 = svdvals(Txx)
ss2 = svdvals(P*Txx)

Plots.plot(title="singular value spectrum")
Plots.plot!(ss1, label="A", l=2)
Plots.plot!(ss2, label="P*A", l=2)

@show norm(u1-u2) / norm(u1)

@show cond(Txx)
@show cond(Tyy)
@show cond(A)

wx, wy, wp = eigvals(Txx), eigvals(Tyy), eigvals(A); println("eigvals found.")

(@isdefined plotresults) || (plotresults = false)
if plotresults
    @eval begin
        using Plots
        plot()
        scatter!(wx,c=:blue,label="primal")
        scatter!(wy,c=:red,label="dual")
        scatter!(wp,c=:green,label="CMP")
        plot!(xlim=(-0.5,0.5),ylim=(-10.0,10.0))
    end
end
