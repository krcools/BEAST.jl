using CompScienceMeshes
using BEAST
using Test

D, Δx = 1.0, 0.35
fn = joinpath(dirname(@__FILE__),"assets","sphere35.in")
#Γ = meshsphere(D, Δx)
Γ = readmesh(fn)
X = raviartthomas(Γ)
Y = subset(X,[1])

o, x, y, z = euclidianbasis(3)
direction, polarisation = z, x
Nt = 200


sol1 = 1.0
Δt1= 0.08/sol1
U1 = timebasisdelta(Δt1, Nt)
W1 = Y ⊗ U1
duration1, delay1, amplitude1 = 8.0/sol1, 12.0/sol1, 1.0
gaussian1 = creategaussian(duration1, delay1, duration1)
E1 = BEAST.planewave(polarisation, direction, derive(gaussian1), sol1)
b1 = assemble(E1,W1)
taxis1 = collect((0:Nt-1)*Δt1)


sol2 = 36.0
Δt2 = 0.08/sol2
U2 = timebasisdelta(Δt2, Nt)
W2 = Y ⊗ U2
duration2, delay2, amplitude2 = 8.0/sol2, 12.0/sol2, 1.0
gaussian2 = creategaussian(duration2, delay2, duration2)
E2 = BEAST.planewave(polarisation, direction, derive(gaussian2), sol2)
b2 = assemble(E2,W2)
taxis2 = collect((0:Nt-1)*Δt2)

gaussian1.(taxis1)[120]
gaussian2.(taxis2)[120]
@test all(gaussian1.(taxis1) .≈ gaussian2.(taxis2))

derive(gaussian1).(taxis1)[120]
derive(gaussian2).(taxis2)[120]/sol2
@test all(derive(gaussian1).(taxis1)/sol1 .≈ derive(gaussian2).(taxis2)/sol2)

timeels1, timead1 = BEAST.assemblydata(U1)
@test timeels1[100][1][1] ≈ Δt1*100
@test timeels1[100][2][1] ≈ Δt1*101

timeels2, timead2 = BEAST.assemblydata(U2)
@test timeels2[100][1][1] ≈ Δt2*100
@test timeels2[100][2][1] ≈ Δt2*101

b1[120]/sol1
b2[120]/sol2
@test all(b1/sol1 .≈ b2/sol2)
