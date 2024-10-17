using CompScienceMeshes
using BEAST

# Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
Γ = meshsphere(radius=1.0, h=0.2)
@show length(Γ)
# X = raviartthomas(Γ)
X = BEAST.gwpdiv(Γ; order=2)

κ, η = 3.0, 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

BEAST.@defaultquadstrat (t,X,X) BEAST.DoubleNumSauterQstrat(2,3,4,4,4,4)
BEAST.@defaultquadstrat (e,X) BEAST.SingleNumQStrat(12)

@hilbertspace j
@hilbertspace k
# efie = @discretise t[k,j]==e[k]  j∈X k∈X
# u, ch = BEAST.gmres_ch(efie; restart=1500)

A = assemble(t[k,j], X, X)
b = assemble(e[k], X)
Ai = BEAST.GMRESSolver(A; reltol=1e-6, restart=1500)
u = Ai * b

include("utils/postproc.jl")
include("utils/plotresults.jl")

# Id = BEAST.Identity()