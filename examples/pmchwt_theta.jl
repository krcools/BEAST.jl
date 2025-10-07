using CompScienceMeshes, BEAST

using LinearAlgebra
using PlotlyJS

p=1
h=0.5
M = meshsphere(1.0, h; generator=:gmsh)
X = BEAST.gwpdiv(M;order=p)
Y = BEAST.gwpdiv(M;order=p+2)

κ,  η  = 1.0, 1.0
κ′, η′ = √2.0κ, η/√2.0
α, α′ = 1/η, 1/η′

T = Maxwell3D.singlelayer(wavenumber=κ);
T′ = Maxwell3D.singlelayer(wavenumber=κ′)
K  = Maxwell3D.doublelayer(wavenumber=κ)
K′ = Maxwell3D.doublelayer(wavenumber=κ′)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ);
H = -1/(im*κ*η)*curl(E)

e = (n × E) × n
h = (n × H) × n

BEAST.@defaultquadstrat (T,X,X) BEAST.DoubleNumSauterQstrat(7,7,6,6,6,6)
BEAST.@defaultquadstrat (K,X,X) BEAST.DoubleNumSauterQstrat(7,7,6,6,6,6)
BEAST.@defaultquadstrat (e,X) BEAST.SingleNumQStrat(10)
BEAST.@defaultquadstrat (h,X) BEAST.SingleNumQStrat(10)

ΘΣ = BEAST.ThetaStars()
ΘΛ = BEAST.ThetaLoops()

@hilbertspace j m
@hilbertspace r s

Θ = ΘΣ[r,j] + ΘΛ[r,j] + ΘΣ[s,m] + ΘΛ[s,m] 

Θh = assemble(Θ,s∈Y, r∈Y, j∈X, m∈X)

Ah = assemble( (η*T+η′*T′)[r,j] -      (K+K′)[r,m] +
                    (K+K′)[s,j] + (α*T+α′*T′)[s,m], s∈X, r∈X, j∈X, m∈X)

Ch = assemble(T[r,j]+T[s,m],s∈Y, r∈Y, j∈Y, m∈Y)

bh = assemble(-e[j]-h[m],j ∈ X, m∈X)

u, stats = solve(BEAST.GMRES(Ah; M=Θh'*Ch*Θh, restart=true, atol=1e-8, rtol=1e-8, verbose=1, memory=50),bh)

uref, statsref = solve(BEAST.GMRES(Ah; restart=true, atol=1e-8, rtol=1e-8, verbose=1, memory=500),bh)

fcrj, _ = facecurrents(u[j],X)
fcrm, _ = facecurrents(u[m],X)

plot(patch(M, norm.(fcrj)),Layout(title="j PMCHWT"))
plot(patch(M, norm.(fcrm)),Layout(title="m PMCHWT"))
