using CompScienceMeshes
using BEAST
using LinearForms

include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))

x = point(1.0, 0.0, 0.0)
y = point(0.0, 1.0, 0.0)
z = point(0.0, 0.0, 1.0)
ι = complex(0.0, 1.0)

Γ = meshsphere(1.0, 0.2)
X = raviartthomas(Γ)

ω = 1.0
μ = ϵ = 1.0

κ = ω * √(ϵ*μ)
η = √(μ/ϵ)

T = MWSingleLayer3D(im*κ)
K = MWDoubleLayer3D(im*κ)

E = PlaneWaveMW(z, x, κ, 1.0)
H = -1/(ι*μ*ω)*curl(E)

e = (n × E) × n
h = (n × H) × n

j, m, = hilbertspace(:j, :m)
k, l, = hilbertspace(:k, :l)

PMCHWT = @varform begin
    2T[k,j] + 2K[k,m] -
    2K[l,j] + 2T[l,m] ==
    h[j] - e[m]
end

pmchwt = @discretise PMCHWT j∈X m∈X k∈X l∈X
u = solve(pmchwt)

u1 = u[1:numfunctions(X)]
fcr, geo = facecurrents(u1,X)
#A = [real(norm(f)) for f in fcr]
p = patch(geo, real.(norm.(fcr)))
PlotlyJS.plot(p)
