

struct Gaussian{T}
    scaling::T
    width::T
    delay::T
end

Gaussian(;scaling=1.0, width, delay) = Gaussian(scaling, width, delay)

(g::Gaussian)(s::Real) = 4*g.scaling/(g.width*√π) * exp(-(4*(s-g.delay)/g.width)^2)


function creategaussian(width,s0,scaling=one(typeof(width)))
    #f(s) = 4*scaling/(width*sqrt(π)) * exp(-(4*(s-s0)/width)^2)
    Gaussian(scaling, width, s0)
    #f(s) = scaling * exp(-(4*(s-s0)/width)^2)
end


function  fouriertransform(g::Gaussian)
    scaling = g.scaling
    width = g.width
    s0 = g.delay
    ft(w) = scaling * exp(-im*w*s0 - (width*w/8)^2) / sqrt(2π)
end


function fouriertransform(a::Array, dt, t0, dim=1)
    n = size(a,dim)
    dω = 2π / (n*dt)
    b = fftshift(fft(a, dim), dim) * dt / sqrt(2π)
    ω0 = -dω * div(n,2)
    b, dω, ω0
end

fouriertransform(a::Array; stepsize, offset, dim=1) = fouriertransform(a, stepsize, offset, dim)


derive(g::Gaussian) =  s -> g(s) * (-8 * (s-g.delay)/g.width) * (4/g.width)


struct ErrorFunction{T}
    scaling::T
    width::T
    delay::T
end

function (f::ErrorFunction)(s)
    f.scaling * 0.5 * (1 + erf(4*(s-f.delay)/f.width))
end

function integrate(f::Gaussian)
    return ErrorFunction(f.scaling, f.width, f.delay)
end
