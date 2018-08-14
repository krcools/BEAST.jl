# If you want a scatter plot of the spectra, define `plotresults = true`
# prior to running this script.

using CompScienceMeshes, BEAST

Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
println("Mesh with $(numvertices(Γ)) vertices and $(numcells(Γ)) cells.")
X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

κ = ω = 1.0; γ = κ*im
T = MWSingleLayer3D(γ)
N = NCross()

Txx = assemble(T,X,X); println("primal discretisation assembled.")
Tyy = assemble(T,Y,Y); println("dual discretisation assembled.")
Nxy = assemble(N,X,Y); println("duality form assembled.")

iNxy = inv(Nxy); println("duality form inverted.")
A = iNxy' * Tyy * iNxy * Txx

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
