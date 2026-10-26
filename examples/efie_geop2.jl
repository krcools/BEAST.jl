using CompScienceMeshes
using BEAST
using PlotlyJS
using LinearAlgebra
const CSM = CompScienceMeshes

κ, η = 3.0, 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

radius = 1.0
h = 0.3

Γ = CSM.gmshsphere(radius, h; order=1)
X = BEAST.gwpdiv(Γ; order=0)

BEAST.@defaultquadstrat (t,X,X) BEAST.DoubleNumSauterQstrat(4,5,6,6,6,6)
BEAST.@defaultquadstrat (e,X) BEAST.SingleNumQStrat(12)

@hilbertspace j
@hilbertspace k

A = assemble(t[k,j], ∏(X), ∏(X); threading=:cellcoloring)
b = assemble(e[k], ∏(X))
Ai = BEAST.lu(AbstractMatrix(A))
u = Ai * b

fcr1, geo1 = facecurrents(u, X; order=3)
tr = patch(geo1, norm.(fcr1); opacity=1.0)
p1 = Plot([tr])

Γ = CSM.gmshsphere(radius, h; order=2)
X = BEAST.gwpdiv(Γ; order=2)

A = assemble(t[k,j], ∏(X), ∏(X); threading=:cellcoloring)
b = assemble(e[k], ∏(X))
Ai = BEAST.lu(AbstractMatrix(A))
u = Ai * b

fcr2, geo2 = facecurrents(u, X; order=3)
tr = patch(geo2, norm.(fcr2); opacity=1.0)
p2 = Plot([tr])

display([p1; p2])