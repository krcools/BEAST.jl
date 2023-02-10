

struct HH3DPlaneWave{T,P}
    direction::P
    wavenumber::T
    amplitude::T
end

function HH3DPlaneWave(direction, wavenumber; amplitude=1)
    w, a = promote(wavenumber, amplitude)
    HH3DPlaneWave(direction, w, a)
end

function (f::HH3DPlaneWave)(r)
    d = f.direction
    k = f.wavenumber
    a = f.amplitude
    a * exp(-im*k*dot(d,r))
end

"""
    HH3DLinearPotential

A potential that linearly increases in `direction` with scaling coefficient `amplitude`.
Its negative gradient will be a uniform vector field pointing in the opposite direction.
"""
struct HH3DLinearPotential{T,P}
    direction::P
    amplitude::T
end

function HH3DLinearPotential(; direction=SVector(1,0,0), amplitude=1.0)
    HH3DLinearPotential(direction ./ norm(direction), amplitude)
end

function (f::HH3DLinearPotential)(r)
    d = f.direction
    a = f.amplitude
    return a * dot(d, r)
end

struct gradHH3DLinearPotential{T,P}
    direction::P
    amplitude::T
end

function gradHH3DLinearPotential(;direction=SVector(0.0,0.0,0.0), amplitude=1.0)
    gradHH3DLinearPotential(direction, amplitude)
end

function (f::gradHH3DLinearPotential)(r)
    d = f.direction
    a = f.amplitude

    return a * d
end

function grad(m::HH3DLinearPotential)
    return gradHH3DLinearPotential(m.direction, m.amplitude)
end

*(a::Number, m::HH3DLinearPotential) = HH3DLinearPotential(m.direction, a * m.amplitude)
*(a::Number, m::gradHH3DLinearPotential) = gradHH3DLinearPotential(m.direction, a * m.amplitude)

dot(::NormalVector, m::gradHH3DLinearPotential) = NormalDerivative(HH3DLinearPotential(m.direction, m.amplitude))

"""
    HH3DMonopole

Potential of a monopole-type point source (e.g., of an electric charge)
"""
struct HH3DMonopole{T,P}
    position::P
    wavenumber::T
    amplitude::T
end

function HH3DMonopole(;position=SVector(0.0,0.0,0.0), wavenumber=0.0, amplitude=1.0)
    w, a = promote(wavenumber, amplitude)
    HH3DMonopole(position, w, a)
end

function (f::HH3DMonopole)(r)
    k = f.wavenumber
    p = f.position
    a = f.amplitude

    return  a*exp(-im * k * norm(r - p)) / (norm(r - p))
end

struct gradHH3DMonopole{T,P}
    position::P
    wavenumber::T
    amplitude::T
end

function gradHH3DMonopole(;position=SVector(0.0,0.0,0.0), wavenumber=0.0, amplitude=1.0)
    w, a = promote(wavenumber, amplitude)
    gradHH3DMonopole(position, w, a)
end

function (f::gradHH3DMonopole)(r)
    a = f.amplitude
    k = f.wavenumber
    p = f.position
    vecR = r - p
    R = norm(vecR)

    return -a * vecR * exp(-im * k * R) / R^2 * (im * k + 1 / R)
end

function grad(m::HH3DMonopole)
    return gradHH3DMonopole(m.position, m.wavenumber, m.amplitude)
end

*(a::Number, m::HH3DMonopole) = HH3DMonopole(m.position, m.wavenumber, a * m.amplitude)
*(a::Number, m::gradHH3DMonopole) = gradHH3DMonopole(m.position, m.wavenumber, a * m.amplitude)

dot(::NormalVector, m::gradHH3DMonopole) = NormalDerivative(HH3DMonopole(m.position, m.wavenumber, m.amplitude))

mutable struct DirichletTrace{F} <: Functional
    field::F
end

function (ϕ::DirichletTrace)(p)
    F = ϕ.field
    x = cartesian(p)
    return F(x)
end

integrand(::DirichletTrace, test_vals, field_vals) = dot(test_vals[1], field_vals)

struct NormalDerivative{F} <: Functional
    field::F
end

const ∂n = Val{:normalderivative}
(::Type{Val{:normalderivative}})(f) = NormalDerivative(f)

function (f::NormalDerivative{T})(manipoint) where T<:HH3DPlaneWave
    d = f.field.direction
    k = f.field.wavenumber
    a = f.field.amplitude
    n = normal(manipoint)
    r = cartesian(manipoint)
    -im*k*a * dot(d,n) * exp(-im*k*dot(d,r))
end

function (f::NormalDerivative{T})(manipoint) where T<:HH3DLinearPotential
    gradient = f.field.amplitude * f.field.direction
    n = normal(manipoint)
    return dot(n, gradient)
end

function (f::NormalDerivative{T})(manipoint) where T<:HH3DMonopole
    m = f.field
    grad_m = grad(m)
    n = normal(manipoint)
    r = cartesian(manipoint)

    return dot(n, grad_m(r))
end

integrand(::NormalDerivative, test_vals, field_vals) = dot(test_vals[1], field_vals)