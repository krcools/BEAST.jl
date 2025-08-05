using CompScienceMeshes
using BEAST

κ1 = 1.0
κ2 = 2.0

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ1)
H = -1/(im*κ1)*curl(E)

T1  = Maxwell3D.singlelayer(wavenumber=κ1)
K1  = Maxwell3D.doublelayer(wavenumber=κ1)
T2  = Maxwell3D.singlelayer(wavenumber=κ2)
K2  = Maxwell3D.doublelayer(wavenumber=κ2)

e = (n × E) × n
h = (n × H) × n

@hilbertspace j m
@hilbertspace k l

a =
    (T1+T2)[k,j] + (-K1-K2)[k,m] +
    (K1+K2)[l,j] + (T1+T2)[l,m]
b =
    -e[k] -h[l]

Γ = meshsphere(;radius=1.0, h=0.35)
RT = raviartthomas(Γ)

X = BEAST.DirectProductSpace([RT,RT])

bx = assemble(b, X)
Axx = assemble(a, X, X)

SXX = BEAST.GMRESSolver(Axx)
uX = SXX * bx
typeof(uX)

import PlotlyBase
import PlotlyDocumenter # hide
using LinearAlgebra

fcrm, geom = facecurrents(uX[m], RT)
fcrj, geoj = facecurrents(uX[j], RT)

ptm = CompScienceMeshes.patch(geom, norm.(fcrm); caxis=(0,1.2) , showscale=false)
ptj= CompScienceMeshes.patch(geoj, norm.(fcrj); caxis=(0,1.2))

pl = [ PlotlyBase.Plot(ptm) PlotlyBase.Plot(ptj) ]

Ty  = Maxwell3D.singlelayer(gamma=κ1)

c = Ty[k,j] + Ty[l,m]

Nx = BEAST.NCross()
d = Nx[k,j] + Nx[l,m]

BC = buffachristiansen(Γ)
Y = BEAST.DirectProductSpace([BC,BC])

Cyy = assemble(c, Y, Y)
Dxy = assemble(d, X, Y)

Di = BEAST.GMRES(Dxy, verbose=0)
Dt = BEAST.GMRES(Dxy', verbose=0)
P = Dt * Cyy * Di

S1 = BEAST.GMRES(Axx; restart=true, atol=1e-8, rtol=1e-8, verbose=0, memory=200)
u1, ch1 = solve(S1, bx);

S2 = BEAST.GMRES(Axx; M=P, restart=true, atol=1e-6, rtol=1e-6, verbose=1, memory=200)
u2, ch2 = solve(S2, bx);

using LinearAlgebra
norm(u1-u2), ch1.niter, ch2.niter
