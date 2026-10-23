using CompScienceMeshes
using BEAST
const CSM = CompScienceMeshes

radius = 1.0
h = 0.2
# Γ = meshsphere(radius=1.0, h=0.2)
Γ = CSM.gmshsphere(radius, h; order=2)
ℰ = CSM.skeleton(Γ, 1)
@assert dimension(Γ) == 2
X = BEAST.gwpdiv(Γ; order=2)
@assert numfunctions(X) == 3 * length(ℰ) + 6 * length(Γ)

κ, η = 3.0, 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

BEAST.@defaultquadstrat (t,X,X) BEAST.DoubleNumSauterQstrat(2,3,4,4,4,4)
BEAST.@defaultquadstrat (e,X) BEAST.SingleNumQStrat(12)

@hilbertspace j
@hilbertspace k

A = assemble(t[k,j], X, X)
b = assemble(e[k], X)
Ai = BEAST.GMRESSolver(A; reltol=1e-6, restart=1500)
u = Ai * b


include("utils/postproc.jl")
include("utils/plotresults.jl")

# Id = BEAST.Identity()