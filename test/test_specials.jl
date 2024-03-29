using BEAST

for T in [Float32, Float64]
    width, delay = T(4.0), T(6.0)
    f = BEAST.Gaussian(width=width, delay=delay)
    g = BEAST.integrate(f)

    step = width/150
    x = range(delay-2*width, stop=delay+2*width, step=step)
    # xc = 0.5*(x[1:end-1] + x[2:end])

    y1 = g.(x)
    y2 = cumsum(f.(x))*step

    using LinearAlgebra
    @assert norm(y1-y2, Inf) < T(1e-2)
end