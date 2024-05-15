struct HH3DPlaneWave{P,K,T}
    direction::P
    gamma::K
    amplitude::T
end

function (f::HH3DPlaneWave)(r)
    d = f.direction
    γ = f.gamma
    a = f.amplitude
    return a * exp(-γ*dot(d,r))
end

function (f::HH3DPlaneWave)(r::CompScienceMeshes.MeshPointNM)
    return f(cartesian(r))
end

scalartype(f::HH3DPlaneWave{P,K,T}) where {P,K,T} = promote_type(eltype(P), K, T)

"""
    HH3DLinearPotential

A potential that linearly increases in `direction` with scaling coefficient `amplitude`.
Its negative gradient will be a uniform vector field pointing in the opposite direction.
"""
struct HH3DLinearPotential{P,T}
    direction::P
    amplitude::T
end

scalartype(f::HH3DLinearPotential{P,T}) where {P,T} = promote_type(eltype(P), T)

function (f::HH3DLinearPotential)(r)
    d = f.direction
    a = f.amplitude
    return a * dot(d, r)
end

function (f::HH3DLinearPotential)(r::CompScienceMeshes.MeshPointNM)
    return f(cartesian(r))
end
struct gradHH3DLinearPotential{T,P}
    direction::P
    amplitude::T
end

function (f::gradHH3DLinearPotential)(r)
    d = f.direction
    a = f.amplitude

    return a * d
end

function (f::gradHH3DLinearPotential)(r::CompScienceMeshes.MeshPointNM)
    r = cartesian(mp)
    return dot(normal(mp), f(r))
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
struct HH3DMonopole{P,K,T}
    position::P
    gamma::K
    amplitude::T
end

scalartype(x::HH3DMonopole{P,K,T}) where {P,K,T} = promote_type(eltype(P), K, T)

function (f::HH3DMonopole)(r::CompScienceMeshes.MeshPointNM)
    return f(cartesian(r))
end

function (f::HH3DMonopole)(r)
    γ = f.gamma
    p = f.position
    a = f.amplitude

    return  a*exp(-γ * norm(r - p)) / (norm(r - p))
end

struct gradHH3DMonopole{P,K,T}
    position::P
    gamma::K
    amplitude::T
end

scalartype(x::gradHH3DMonopole{P,K,T}) where {P,K,T} = promote_type(eltype(P), K, T)

function (f::gradHH3DMonopole)(r)
    a = f.amplitude
    γ = f.gamma
    p = f.position
    vecR = r - p
    R = norm(vecR)

    return -a * vecR * exp(-γ * R) / R^2 * (γ + 1 / R)
end

function (f::gradHH3DMonopole)(mp::CompScienceMeshes.MeshPointNM)
    r = cartesian(mp)
    return dot(normal(mp), f(r))
end

function grad(m::HH3DMonopole)
    return gradHH3DMonopole(m.position, m.gamma, m.amplitude)
end

*(a::Number, m::HH3DMonopole) = HH3DMonopole(m.position, m.gamma, a * m.amplitude)
*(a::Number, m::gradHH3DMonopole) = gradHH3DMonopole(m.position, m.gamma, a * m.amplitude)

dot(::NormalVector, m::gradHH3DMonopole) = NormalDerivative(HH3DMonopole(m.position, m.gamma, m.amplitude))

mutable struct DirichletTrace{T,F} <: Functional
    field::F
end

DirichletTrace(f::F) where {F} = DirichletTrace{scalartype(f), F}(f)
DirichletTrace{T}(f::F) where {T,F} = DirichletTrace{T,F}(f)
scalartype(s::DirichletTrace{T}) where {T} = T

function (ϕ::DirichletTrace)(p)
    F = ϕ.field
    x = cartesian(p)
    return F(x)
end

integrand(::DirichletTrace, test_vals, field_vals) = dot(test_vals[1], field_vals)

struct NormalDerivative{T,F} <: Functional
    field::F
end

NormalDerivative(f::F) where {F} = NormalDerivative{scalartype(f), F}(f)
NormalDerivative{T}(f::F) where {T,F} = NormalDerivative{T,F}(f)
scalartype(s::NormalDerivative{T}) where {T} = T

const ∂n = Val{:normalderivative}
(::Type{Val{:normalderivative}})(f) = NormalDerivative(f)

function (f::NormalDerivative)(nbd)
    x = cartesian(nbd)
    g = gradient(f.field)(x)
    n = normal(nbd)
    return dot(n,g)
end

function (f::NormalDerivative{T,F})(manipoint) where {T,F<:HH3DPlaneWave}
    d = f.field.direction
    γ = f.field.gamma
    a = f.field.amplitude
    n = normal(manipoint)
    r = cartesian(manipoint)
    return -γ*a * dot(d,n) * exp(-γ*dot(d,r))
end

function (f::NormalDerivative{T,F})(manipoint) where {T,F<:HH3DLinearPotential}
    gradient = f.field.amplitude * f.field.direction
    n = normal(manipoint)
    return dot(n, gradient)
end

function (f::NormalDerivative{T,F})(manipoint) where {T,F<:HH3DMonopole}
    m = f.field
    grad_m = grad(m)
    n = normal(manipoint)
    r = cartesian(manipoint)

    return dot(n, grad_m(r))
end

integrand(::NormalDerivative, test_vals, field_vals) = dot(test_vals[1], field_vals)