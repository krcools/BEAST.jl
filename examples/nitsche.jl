using CompScienceMeshes, BEAST
o, x, y, z = euclidianbasis(3)

trc = X->ntrace(X,γ)
dvg = divergence

h1 = h2 = h3 = 1/15; h = max(h1, h2, h3)
width, height = 1.0, 0.5
Γ1 = meshrectangle(width, height, h1)
Γ2 = meshrectangle(width, height, h2); CompScienceMeshes.rotate!(Γ2, 0.5π * [1,0,0])
Γ3 = meshrectangle(width, height, h3); CompScienceMeshes.rotate!(Γ3, 1.0π * [1,0,0])
γ = meshsegment(width, width, 3)

κ = 1.0
S, T, I = SingleLayerTrace(κ), MWSingleLayer3D(κ), Identity()
E = planewavemw3d(direction=z, polarization=x, wavenumber=κ)
e = (n×E)×n

X1, X2, X3 = raviartthomas(Γ1, γ), raviartthomas(Γ2, γ), raviartthomas(Γ3, γ)
X = X1 × X2 × X3

@hilbertspace j
@hilbertspace k
α, β = 1/(im*κ), log(abs(h))/(im*κ)
Eq = @varform T[k,j] + α*S[trc(k), dvg(j)] + α*S.'[dvg(k), trc(j)] + β*I[trc(k), trc(j)] == e[k]
eq = @discretise Eq  j∈X k∈X
u = solve(eq)

fcr, geo = facecurrents(u, X)
include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))
p = patch(geo, real.(norm.(fcr)))
PlotlyJS.plot(p)
