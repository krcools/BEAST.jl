using CompScienceMeshes, BEAST
using Test

m = meshsphere(1.0,0.8)
X = BEAST.raviartthomas(m)
c = 1.0 #speed of light
Δt, Nt = 0.2, 10
V = BEAST.FiniteDiffTimeStep(X, Δt, Nt)
G = BEAST.assemble(Identity(), X, X)

LaplaceId(s::T) where {T} = s*Identity()
kmax = 5
rho = 1.0001


method = BEAST.BE(Δt)
T = BEAST.FiniteDiffConvolutionQuadrature(LaplaceId, method, Δt, kmax, rho)
Z = BEAST.assemble(T, V, V)
@test Z.data[:,:,1]*Δt ≈ G
@test Z.data[:,:,2]*Δt*(-1) ≈ G
for k in 3:kmax
    @test norm(Z.data[:,:,k]) ≈ 0 atol=1e-10 
end

method = BEAST.BDF2(Δt)
T = BEAST.FiniteDiffConvolutionQuadrature(LaplaceId, method, Δt, kmax, rho)
Z = BEAST.assemble(T, V, V)
@test Z.data[:,:,1]*Δt*2/3 ≈ G
@test Z.data[:,:,2]*Δt*(-2)/4 ≈ G
@test Z.data[:,:,3]*Δt*(2)/1 ≈ G
for k in 4:kmax
    @test norm(Z.data[:,:,k]) ≈ 0 atol=1e-10 
end