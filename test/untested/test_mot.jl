using BEAST
using Test

Id = ones(1,1)
Δt = 0.05
Z = zeros(1,1,2)
Z[:,:,1] = Id
Z[:,:,2] = Δt - Id
#Z = [Id; (Δt-Id)]

Nt = 100
B = zeros(1,Nt)
B[:,1] = ones(1)

Z0 = Z[:,:,1]
W0 = inv(Z0)

x = BEAST.mot(W0,Z,B,Nt)

n = 20
@show x[n]
@show exp(-n*Δt)
@show norm(x[n]-exp(-n*Δt))

@test norm(x[n]-exp(-n*Δt)) < 1.0e-2
