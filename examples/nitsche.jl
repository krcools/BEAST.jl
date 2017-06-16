using Base.Test
using LinearForms
using CompScienceMeshes
using BoundaryElements

x = point(1.0, 0.0, 0.0)
y = point(0.0, 1.0, 0.0)
z = point(0.0, 0.0, 1.0)

n = BoundaryElements.n

# Define the kernels
κ = 1.0
S = SingleLayerTrace(κ)
T = MWSingleLayer3D(κ)
I = Identity()

# Build the geometry
h1, h2 = 1/12, 1/15
h = max(h1,h2)
W, H = 1.0, 0.5
Γ1 = meshrectangle(W,H,h1)
Γ2 = translate(meshrectangle(W,H,h2), point(0,-H,0))
γ = meshsegment(W,W,3)

X1 = raviartthomas(Γ1, γ)
X2 = raviartthomas(Γ2, γ)
X = X1 × X2

# Define the incident field
E = PlaneWaveMW(z, y, κ, complex(1.0))
e = (n×E)×n

j, = hilbertspace(:j)
k, = hilbertspace(:k)

trc = X->ntrace(X,γ)
dvg = divergence

α = 1/(im*κ)
β = α*log(abs(h))

Eq = @varform begin
    T[k,j]              +
    α*S[trc(k), dvg(j)]   +
    α*S.'[dvg(k), trc(j)] +
    β*I[trc(k), trc(j)] == e[k]
end

eq = @discretise Eq j∈X k∈X

b = rhs(eq)
Z = sysmatrix(eq)

u = Z \ b

## Post-processing
fcr, geo = facecurrents(u, X)

include(Pkg.dir("CompScienceMeshes","examples","matlab_patches.jl"))
p = patch(geo, real.(norm.(fcr)))
PlotlyJS.plot(p)
