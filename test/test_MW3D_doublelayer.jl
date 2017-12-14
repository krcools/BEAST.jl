using BEAST
using CompScienceMeshes

Γ = readmesh(Pkg.dir("BEAST", "test","assets","flatplate.in"))
#Γ = readmesh(Pkg.dir("BEAST", "test","assets","smallflatplate.in"))
RT = raviartthomas(Γ)

X = RT
Y = RT

κ = 1.0
t = BEAST.MWDoubleLayer3D(κ*im)

BEAST.quadrule(op::BEAST.MaxwellOperator3D, g::BEAST.RTRefSpace, f::BEAST.RTRefSpace, i, τ, j, σ, qd) = BEAST.q(op, g, f, i, τ, j, σ, qd)
Zss = assemble(t,X,Y)

BEAST.quadrule(op::BEAST.MaxwellOperator3D, g::BEAST.RTRefSpace, f::BEAST.RTRefSpace, i, τ, j, σ, qd) = BEAST.qrdf(op, g, f, i, τ, j, σ, qd)
Zqrdf = assemble(t,X,Y)

BEAST.quadrule(op::BEAST.MaxwellOperator3D, g::BEAST.RTRefSpace, f::BEAST.RTRefSpace, i, τ, j, σ, qd) = BEAST.qrib(op, g, f, i, τ, j, σ, qd)
Zbogaert = assemble(t,X,Y)

println(abs(Zss)./abs(Zbogaert))
