using CompScienceMeshes
using BEAST
using Test
import BEAST: BasisFunction, ∇ 

m = meshrectangle(1.0, 1.0, 0.5, 3)
k = 1.2

K = Maxwell3D.doublelayer(wavenumber = k)
I = Identity()
r = raviartthomas(m)

KM = assemble(K,r,r)
IM = assemble(I,r,r)
CM = assemble(NCross(),r,r) # normal on test function
#limit from inside
TM = KM +1/2*CM


G = BEAST.greenhh3d(wavenumber=k)
KC = BEAST.build_potential(∇(G)×BasisFunction())
KCT1 = BEAST.γₜ(KC,-1)
TCM1 = assemble(KCT1,r,r)

@test TM ≈ TCM1 atol=sqrt(eps())

#limit from outside
TM = KM -1/2*CM


G = BEAST.greenhh3d(wavenumber=k)
KC = BEAST.build_potential(∇(G)×BasisFunction(),m)
KCT1 = BEAST.γₜ(KC,1)
TCM1 = assemble(KCT1,r,r)


@test TM ≈ TCM1 atol=sqrt(eps())

## same for non matching meshes

m1 = meshrectangle(1.0, 1.0, 0.5, 3)
m2 = translate(m1,[0.2,0.2,0.0])
k = 1.2

K = Maxwell3D.doublelayer(wavenumber = k)
I = Identity()
r1 = raviartthomas(m1)
r2 = raviartthomas(m2)

KM = assemble(K,r1,r2)
IM = assemble(I,r1,r2)
CM = assemble(NCross(),r1,r2) # normal on test function
#limit from inside
TM = KM +1/2*CM


G = BEAST.greenhh3d(wavenumber=k)
KC = BEAST.build_potential(∇(G)×BasisFunction(),m2)
KCT1 = BEAST.γₜ(KC,-1)
TCM1 = assemble(KCT1,r1,r2)


@test TM ≈ TCM1 atol=sqrt(eps())

#limit from outside
TM = KM -1/2*CM


G = BEAST.greenhh3d(wavenumber=k)
KC = BEAST.build_potential(∇(G)×BasisFunction(),m2)
KCT1 = BEAST.γₜ(KC,1)
TCM1 = assemble(KCT1,r1,r2)


@test TM ≈ TCM1 atol=sqrt(eps())

