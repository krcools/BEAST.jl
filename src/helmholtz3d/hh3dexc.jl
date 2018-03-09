

struct HH3DPlaneWave{T,P} <: Functional
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


integrand(::NormalDerivative, test_vals, field_vals) = dot(test_vals[1], field_vals)
