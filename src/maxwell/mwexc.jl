mutable struct PlaneWaveMW{T,P}
  direction::P
  polarisation::P
  wavenumber::T
  amplitude::T
end

function PlaneWaveMW(d,p,k,a = 1)
    T = promote_type(eltype(d), eltype(p), typeof(k), typeof(a))
    P = similar_type(typeof(d), T)
    PlaneWaveMW{T,P}(d,p,k,a)
end

"""
    planewavemw3d(;direction, polarization, wavenumber[, amplitude=1])

Create a plane wave solution to Maxwell's equations.
"""
planewavemw3d(;
    direction    = error("missing arguement `direction`"),
    polarization = error("missing arguement `polarization`"),
    wavenumber   = error("missing arguement `wavenumber`"),
    amplitude    = 1,
    ) = PlaneWaveMW(direction, polarization, wavenumber, amplitude)

function (e::PlaneWaveMW)(x)
  k = e.wavenumber
  d = e.direction
  u = e.polarisation
  a = e.amplitude
  a * exp(-im * k * dot(d, x)) * u
end

function curl(field::PlaneWaveMW)
  k = field.wavenumber
  d = field.direction
  u = field.polarisation
  a = field.amplitude
  v = d × u
  b = -a * im * k
  PlaneWaveMW(d, v, k, b)
end

*(a::Number, e::PlaneWaveMW) = PlaneWaveMW(e.direction, e.polarisation, e.wavenumber, a*e.amplitude)

mutable struct CrossTraceMW{F} <: Functional
  field::F
end

mutable struct TangTraceMW{F} <: Functional
  field::F
end

cross(::NormalVector, p::Function) = CrossTraceMW(p)
cross(t::CrossTraceMW, ::NormalVector) = TangTraceMW(t.field)

function (ϕ::CrossTraceMW)(p)
  F = ϕ.field
  x = cartesian(p)
  n = normal(p)
  return n × F(x)
end

function (ϕ::TangTraceMW)(p)
  F = ϕ.field
  x = cartesian(p)
  n = normal(p)
  return (n × F(x)) × n
end

integrand(::TangTraceMW, gx, ϕx) = gx[1] ⋅ ϕx
integrand(::CrossTraceMW, test_vals, field_val) = test_vals[1] ⋅ field_val

struct NDotTrace{F} <: Functional
  field::F
end

(ϕ::NDotTrace)(p) = dot(normal(p), ϕ.field(cartesian(p)))
integrand(::NDotTrace, g, ϕ) = dot(g.value, ϕ)
LinearAlgebra.dot(::NormalVector, f) = NDotTrace(f)
