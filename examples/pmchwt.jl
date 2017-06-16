using CompScienceMeshes, BEAST

Γ = meshcuboid(1.0, 1.0, 1.0, 0.15)
X = raviartthomas(Γ)

κ,  η  = 1.0, 1.0
κ′, η′ = 1.4κ, η/1.4

γ, γ′ = im*κ, im*κ′
T,  K  = MWSingleLayer3D(γ),  MWDoubleLayer3D(γ)
T′, K′ = MWSingleLayer3D(γ′), MWDoubleLayer3D(γ′)

o, x, y, z = euclidianbasis(3)
d = normalize(x+y+z); p = normalize(d × z)
E = planewavemw3d(direction=d, polarization=p, wavenumber=κ)
H = -1/(im*κ*η)*curl(E)

e, h = (n × E) × n, (n × H) × n

@hilbertspace j m
@hilbertspace k l

α, α′ = 1/η, 1/η′
pmchwt = @discretise(
    (η*T+η′*T′)[k,j] +      (K+K′)[k,m] -
         (K+K′)[l,j] + (α*T+α′*T′)[l,m] == h[k] + e[l],
    j∈X, m∈X, k∈X, l∈X)

u = solve(pmchwt)

## Post-processing
u1 = u[1:numfunctions(X)]
fcr, geo = facecurrents(u1,X)

include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))
p = patch(geo, real.(norm.(fcr)))
PlotlyJS.plot(p)
