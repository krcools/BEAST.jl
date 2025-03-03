using CompScienceMeshes
using BEAST

# Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
# Γ = meshsphere(radius=1.0, h=0.1)
Γ = CompScienceMeshes.meshrectangle_unit_tri_2graded(0.3, 3)
@show length(Γ)
X = raviartthomas(Γ)

κ, η = 6.0, 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

@hilbertspace j
@hilbertspace k

qs = BEAST.DoubleNumWiltonSauterQStrat{Int64, Int64}(10, 11, 10, 11, 8, 8, 8, 8)

U = BEAST.DirectProductSpace([X])
A = assemble(t[k,j], U, U; quadstrat=qs)
b = assemble(e[k], U)

u = AbstractMatrix(A) \ b 
fcr, geo = facecurrents(u, X)

import Plotly
using LinearAlgebra
Plotly.plot(patch(geo, log10.(norm.(fcr))))

efie = @discretise t[k,j]==e[k]  j∈X k∈X
u, ch = BEAST.gmres_ch(efie; restart=1500)

include("utils/postproc.jl")
include("utils/plotresults.jl")

