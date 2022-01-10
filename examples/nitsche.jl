using CompScienceMeshes, BEAST

trc = X->ntrace(X,γ)
dvg = divergence

h1 = h2 = h3 = 1/15; h = max(h1, h2, h3)
width, height = 1.0, 0.5
Γ1 = meshrectangle(width, height, h1)
Γ2 = meshrectangle(width, height, h2); CompScienceMeshes.rotate!(Γ2, 0.5π*x̂)
Γ3 = meshrectangle(width, height, h3); CompScienceMeshes.rotate!(Γ3, 1.0π*x̂)
γ = meshsegment(width, width, 3)

κ = 1.0
S, T, I = SingleLayerTrace(κ*im), MWSingleLayer3D(κ), Identity()
St = transpose(S)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n×E)×n

in_interior1 = CompScienceMeshes.interior_tpredicate(Γ1)
in_interior2 = CompScienceMeshes.interior_tpredicate(Γ2)
in_interior3 = CompScienceMeshes.interior_tpredicate(Γ3)

on_junction = CompScienceMeshes.overlap_gpredicate(γ)
edges1 = skeleton(e -> in_interior1(e) || on_junction(chart(Γ1,e)), Γ1,1)
edges2 = skeleton(e -> in_interior2(e) || on_junction(chart(Γ2,e)), Γ2,1)
edges3 = skeleton(e -> in_interior3(e) || on_junction(chart(Γ3,e)), Γ3,1)

X1 = raviartthomas(Γ1, edges1)
X2 = raviartthomas(Γ2, edges2)
X3 = raviartthomas(Γ3, edges3)
# X1, X2, X3 = raviartthomas(Γ1, γ), raviartthomas(Γ2, γ), raviartthomas(Γ3, γ)
X = X1 × X2 × X3

@hilbertspace j
@hilbertspace k
α, β = 1/(im*κ), log(abs(h))/(im*κ)
Eq = @varform T[k,j] + α*S[trc(k), dvg(j)] + α*St[dvg(k), trc(j)] + β*I[trc(k), trc(j)] == e[k]
eq = @discretise Eq  j∈X k∈X
u = solve(eq)

Θ, Φ = range(0.0,stop=π,length=100), 0.0
ffpoints = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]
farfield = potential(MWFarField3D(κ*im), ffpoints, u, X)

using Plots
using LinearAlgebra
plot(Θ,norm.(farfield))
