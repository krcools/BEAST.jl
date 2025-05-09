using CompScienceMeshes
using BEAST

# Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
# Γ = meshsphere(radius=1.0, h=0.1)
Γ = CompScienceMeshes.meshrectangle_unit_tri_2graded(0.25, 3)
@show length(Γ)
X = raviartthomas(Γ)

κ, η = 4.0, 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ŷ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

@hilbertspace j
@hilbertspace k

qs = BEAST.DoubleNumWiltonSauterQStrat{Int64, Int64}(11, 11, 11, 12, 8, 8, 8, 8)
qsl = BEAST.SingleNumQStrat(13)

U = BEAST.DirectProductSpace([X])
A = assemble(t[k,j], U, U; quadstrat=qs)
b = assemble(e[k], U; quadstrat=qsl)

u = AbstractMatrix(A) \ b 
fcr, geo = facecurrents(u, X)

import Plotly
using LinearAlgebra
Plotly.plot(patch(geo, real.(getindex.(fcr,1)), caxis=(-10,10)))
Plotly.plot(patch(geo, real.(getindex.(fcr,2)), caxis=(-10,10)))
Plotly.plot(patch(geo, norm.(fcr), caxis=(0,10)))

# efie = @discretise t[k,j]==e[k]  j∈X k∈X
# u, ch = BEAST.gmres_ch(efie; restart=1500)

# include("utils/postproc.jl")
# include("utils/plotresults.jl")
Id = Identity()
Nuu = assemble(Id[k,j], U, U)
v = AbstractMatrix(Nuu) \ b
fcr, geo = facecurrents(v, X)

Plotly.plot(patch(geo, real(getindex.(fcr,1))))