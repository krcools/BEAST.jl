using BEAST
using Test
using LinearAlgebra

for T in [Float32, Float64]
    width = T(4.0)
    delay = T(6.0)
    g = BEAST.creategaussian(width, delay)
    G = BEAST.fouriertransform(g)

    dx = width / 15
    x0 = T(0.0)
    x = x0 : dx : 3*delay
    n = length(x)
    y = g.(x)

    Y, dω, ω0 = BEAST.fouriertransform(y, dx, x0)
    ω = collect(ω0 .+ (0:n-1)*dω)
    Z = G.(ω)

    # using Plots
    # plot(ω, real(Y))
    # scatter!(ω, real(Z))
    # plot!(ω, imag(Y))
    # scatter!(ω, imag(Z))
    @test norm(Y-Z) < sqrt(eps(T))
end