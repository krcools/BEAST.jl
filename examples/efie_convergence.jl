using CompScienceMeshes
using BEAST

using BakerStreet
using DataFrames

function payload(;h)
    # Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
    @show h
    Γ = meshsphere(radius=1.0, h=h)
    # @show length(Γ)
    X = raviartthomas(Γ)

    κ, η = 1.0, 1.0
    t = Maxwell3D.singlelayer(wavenumber=κ)
    s = -Maxwell3D.singlelayer(gamma=κ)
    E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
    e = (n × E) × n

    @hilbertspace j
    @hilbertspace k

    A = assemble(t[k,j], X, X)
    S = assemble(s[k,j], X, X)
    b = assemble(e[k], X)

    Ai = BEAST.GMRESSolver(A; restart=1500, abstol=1e-5, maxiter=10_000)
    u = Ai * b
    u_Snorm = real(sqrt(dot(u, S*u)))
    return (;u, X, u_Snorm)
end

α = 0.8
h = collect(0.4 * α.^(0:10))  
# h = collect(logrange(0.4, 0.05, length=8))
@runsims payload h

df = @collect_results

using LinearAlgebra
using Plots
plot(df.h, real.(df.u_Snorm), marker=:.)

# α = df.h[1] / df.h[2]
W1 = reverse(df.u_Snorm)
W2 = Iterators.drop(W1,1)
W3 = Iterators.drop(W1,2)
p = [log((w2-w1)/(w3-w2)) / log(1/α) for (w1,w2,w3) in zip(W1,W2,W3)]

# efie = @discretise t[k,j]==e[k]  j∈X k∈X
# u, ch = BEAST.gmres_ch(efie; restart=1500)

# include("utils/postproc.jl")
# include("utils/plotresults.jl")

