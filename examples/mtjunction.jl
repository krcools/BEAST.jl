using CompScienceMeshes, BEAST

width, height, h = 1.0, 0.5, 0.05
G1 = meshrectangle(width, height, h)
G2 = CompScienceMeshes.rotate(G1, 0.5π * x̂)
G3 = CompScienceMeshes.rotate(G1, 1.0π * x̂)

G12 = weld(G1,-G2)
G23 = weld(G2,-G3)
G31 = weld(G3,-G1)

X12 = raviartthomas(G12)
X23 = raviartthomas(G23)
X31 = raviartthomas(G31)

Y12 = buffachristiansen(G12)
Y23 = buffachristiansen(G23)
Y31 = buffachristiansen(G31)

κ = 1.0
SL = Maxwell3D.singlelayer(wavenumber=κ)
N = NCross()

X = X12 × X23 × X31
Y = Y12 × Y23 × Y31

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

@hilbertspace j
@hilbertspace k
mtefie = @discretise(
    SL[k,j] == e[k],
    j ∈ X, k ∈ X)

# u = gmres(mtefie)
#
# offset = 1; u12 = u[offset:offset+numfunctions(X12)-1]
# offset += numfunctions(X12); u23 = u[offset:offset+numfunctions(X23)-1]
# offset += numfunctions(X23); u31 = u[offset:offset+numfunctions(X31)-1]
#
# fcr12, _ = facecurrents(u12, X12)
# fcr23, _ = facecurrents(u23, X23)
# fcr31, _ = facecurrents(u31, X31)

import PlotlyJS
using LinearAlgebra
# p1 = PlotlyJS.Plot(patch(G12, norm.(fcr12)))
# p2 = PlotlyJS.Plot(patch(G23, norm.(fcr23)))
# p3 = PlotlyJS.Plot(patch(G31, norm.(fcr31)))

PlotlyJS.plot([p1, p2, p3])

Sxx = BEAST.sysmatrix(mtefie)
duality = @discretise N[k,j]==e[k] k∈X j∈Y
Nxy = BEAST.sysmatrix(duality)

bcefie = @discretise SL[k,j] == e[k] j∈Y k∈Y
Syy = BEAST.sysmatrix(bcefie)

iNxy = inv(Nxy)

gmres()
