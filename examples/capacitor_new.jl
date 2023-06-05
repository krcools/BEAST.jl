using BEAST
using CompScienceMeshes
using LinearAlgebra
using SparseArrays

import Plotly
import Plots
Plots.plotly()

function Base.:+(a::S, b::S) where {S<:BEAST.Space}
    @assert geometry(a) == geometry(b)
    S(geometry(a), [a.fns; b.fns], [a.pos; b.pos])
end

l = w = 1.0 #Length and width of capacitor plates
d = 0.1     #seperation of plates
h = 1/60    #size of meshes

Γ₀ = meshrectangle(l,w,h)
Γ₁ = CompScienceMeshes.translate(Γ₀, point(0.0,0.0,d))

γ₀ = meshsegment(l,h,3)
γ₁ = CompScienceMeshes.translate(γ₀, point(0.0,0.0,d))

Γ = weld(Γ₁, Γ₀)

κ = 2pi / 100
V₀ = 1.0
f = ScalarTrace{typeof(V₀)}(p -> V₀)

edges_all = skeleton(Γ, 1)
edges_int = submesh(!in(boundary(Γ)), edges_all)

on_γ₀ = overlap_gpredicate(γ₀)
edges_pt0 = submesh(edges_all) do m, edge
    ch = chart(m, edge)
    return on_γ₀(ch)
end

on_γ₁ = overlap_gpredicate(γ₁)
edges_pt1 = submesh(edges_all) do m, edge
    ch = chart(m, edge)
    return on_γ₁(ch)
end

RT_int = raviartthomas(Γ, edges_int)
RT_pt0 = raviartthomas(Γ, edges_pt0)
RT_pt1 = raviartthomas(Γ, edges_pt1)

verts_pt0_int = submesh(!in(boundary(edges_pt0)), skeleton(edges_pt0,0))
verts_pt1_int = submesh(!in(boundary(edges_pt1)), skeleton(edges_pt1,0))

C0 = connectivity(verts_pt0_int, edges_pt0)
C1 = connectivity(verts_pt1_int, edges_pt1)

X = RT_int
Y0 = RT_pt0 * C0
Y1 = RT_pt1 * C1

# This is ugly and will fail if the two ports are not equal...
D = sparse(reshape([fill(+1.0, length(edges_pt0)); fill(-1.0, length(edges_pt1))],length(edges_pt0)+length(edges_pt1),1))
RT_pt = RT_pt0 + RT_pt1
Z = RT_pt * D

@hilbertspace x y0 y1 z
@hilbertspace ξ η0 η1 ζ

T = Maxwell3D.singlelayer(wavenumber=κ)
trc = X->ntrace(X,γ₁)

efie = @discretise( @varform(
    T[ξ,x]  + T[ξ,y0]  + T[ξ,y1]  + T[ξ,z]  +
    T[η0,x] + T[η0,y0] + T[η0,y1] + T[η0,z] +
    T[η1,x] + T[η1,y0] + T[η1,y1] + T[η1,z] +
    T[ζ,x]  + T[ζ,y0]  + T[ζ,y1]  + T[ζ,z] == f[trc(ζ)]),
    ξ∈X, η0∈Y0, η1∈Y1, ζ∈Z,
    x∈X, y0∈Y0, y1∈Y1, z∈Z)
u = solve(efie)

# assemble(efie.equation.rhs, efie.test_space_dict)
# @enter assemble(efie.equation.rhs, X + Y0 + Y1 + Z)
# assemble(efie.equation.lhs, efie.test_space_dict, efie.trial_space_dict)

# You can access the current coefficients pertaining to subspaces
# using the Hilbert space 'placeholder'
@show length(u[x])
@show length(u[y0])
@show length(u[y1])
@show length(u[z])

S = (((X + Y0) + Y1) + Z)
fcr, geo = facecurrents(u, S)
Plotly.plot(patch(geo, norm.(fcr))) |> display

# Compute the Scalar Potential across the ports. Note how a voltage
# gap of 1V comes out as a 'natural' condition on the solution.
zs = range(-0.5, stop=0.5, length=100)
pts = [point(0.5,0.5,z) for z ∈ zs]
Φ = potential(MWSingleLayerPotential3D(κ), pts, u, S; type=ComplexF64)
Plots.plot(zs, real(Φ), xlabel="height",ylabel="scalar potential",label=false) |> display

# Plot current along γ₀ and γ₁. Note: in this symmetric situation, we 
# expect these to be opposite on corresponding points of the two ports.
# In general this is not true; in fact, the two ports could have different shape and size
traceS0 = ntrace(S, γ₀)
traceS1 = ntrace(S, γ₁)
fcr0, geo0 = facecurrents(u, traceS0)
fcr1, geo1 = facecurrents(u, traceS1)
Plots.plot()
Plots.plot!(imag.(fcr0))
Plots.plot!(imag.(fcr1)) |> display