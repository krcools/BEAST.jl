using CompScienceMeshes
using BEAST

m = meshsphere(radius=1.0, h=0.35)
X = raviartthomas(m)

@hilbertspace m j
@hilbertspace k l

@hilbertspace u
@hilbertspace v

κ = 3.0
η = 1.0
T = Maxwell3D.singlelayer(wavenumber=κ)
K = Maxwell3D.doublelayer(wavenumber=κ)
H = 0.5 * BEAST.NCross()

A = K[k,m] - T[k,j] + T[l,m] + K[l,j]
B = A[u,v]

Bh = @discretise B u∈X×X v∈X×X
Z = assemble(Bh)