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

abstract type Dipole end

mutable struct DipoleMW{T,P} <: Dipole
    location::P
    orientation::P
    wavenumber::T
end

function DipoleMW(l,o,k)
    T = promote_type(eltype(l), eltype(o), typeof(k))
    P = similar_type(typeof(l), T)
    DipoleMW{T,P}(l,o,k)
end

mutable struct curlDipoleMW{T,P} <: Dipole
    location::P
    orientation::P
    wavenumber::T
end

function curlDipoleMW(l,o,k)
    T = promote_type(eltype(l), eltype(o), typeof(k))
    P = similar_type(typeof(l), T)
    curlDipoleMW{T,P}(l,o,k)
end

"""
    dipolemw3d(;location, orientation, wavenumber)

Create an electric dipole solution to Maxwell's equations representing the electric
field part. Implementation is based on (9.18) of Jackson's “Classical electrodynamics”,
with the notable difference that the ``\exp(ikr)`` is used.
"""
dipolemw3d(;
    location    = error("missing arguement `location`"),
    orientation = error("missing arguement `orientation`"),
    wavenumber   = error("missing arguement `wavenumber`"),
    ) = DipoleMW(location, orientation, wavenumber)

function (d::DipoleMW)(x; isfarfield=false)
    k = d.wavenumber
    x_0 = d.location
    p = d.orientation
    if isfarfield
      # postfactor (4*π*im)/k to be consistent with BEAST far field computation
      # and, of course, adapted phase factor exp(im*k*dot(n,x_0)) with
      # respect to (9.19) of Jackson's Classical Electrodynamics
      r = norm(x)
      n = x/r
      return cross(cross(n,1/(4*π)*(k^2*cross(cross(n,p),n))*exp(im*k*dot(n,x_0))*(4*π*im)/k),n)
    else
      r = norm(x-x_0)
      n = (x - x_0)/r
      return 1/(4*π)*exp(-im*k*r)*(k^2/r*cross(cross(n,p),n) + 
              (1/r^3 + im*k/r^2)*(3*n*dot(n,p) - p))
    end
end

function (d::curlDipoleMW)(x; isfarfield=false)
    k = d.wavenumber
    x_0 = d.location
    p = d.orientation
    if isfarfield
      # postfactor (4*π*im)/k to be consistent with BEAST far field computation
      r = norm(x)
      n = x/r
      return (-im*k)*k^2/(4*π)*cross(n,p)*(4*π*im)/k*exp(im*k*dot(n,x_0))
    else
      r = norm(x-x_0)
      n =  (x - x_0)/r
      return -im*(k^3)/(4*π)*cross(n,p)*exp(-im*k*r)/r*(1 + 1/(im*k*r))
    end
end

function curl(d::DipoleMW)
    return curlDipoleMW(d.location, d.orientation, d.wavenumber)
end

*(a::Number, d::DipoleMW) = DipoleMW(d.location, a .* d.orientation, d.wavenumber)
*(a::Number, d::curlDipoleMW) = curlDipoleMW(d.location, a .* d.orientation, d.wavenumber)

mutable struct CrossTraceMW{F} <: Functional
  field::F
end

mutable struct TangTraceMW{F} <: Functional
  field::F
end

cross(::NormalVector, p::Function) = CrossTraceMW(p)
cross(::NormalVector, p::PlaneWaveMW) = CrossTraceMW(p)
cross(::NormalVector, p::Dipole) = CrossTraceMW(p)
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



mutable struct CurlGreen{T,U,V}
  wavenumber::T
  source::U
  position::V
end

function (f::CurlGreen)(x)
  γ = im * f.wavenumber
  r = x-f.position
  R = norm(r)
  g = exp(-γ*R)/(4π*R)
  j = f.source
  return -γ/R * (1 + 1/(γ*R)) * g * (r × j)
end


mutable struct CurlCurlGreen{T,U,V}
  wavenumber::T
  source::U
  position::V
end

cross(::NormalVector, p::CurlGreen) = CrossTraceMW(p)
cross(::NormalVector, p::CurlCurlGreen) = CrossTraceMW(p)

function (f::CurlCurlGreen)(x)
  γ = im * f.wavenumber
  r = x - f.position
  R = norm(r)
  g = exp(-γ*R)/(4π*R)
  j = f.source
  return (-γ^2/R^2 * (r × j) × r + (1/R^2 + γ/R)/R^2 * (3r * dot(r,j) - R^2 * j)) * g
end

curl(f::CurlGreen) = CurlCurlGreen(f.wavenumber, f.source, f.position)
function curl(f::CurlCurlGreen)
  κ = f.wavenumber
  j = κ^2 * f.source
  x = f.position
  return CurlGreen(κ, j, x)
end

Base.:*(a::Number, f::CurlGreen) = CurlGreen(f.wavenumber, a*f.source, f.position)
Base.:*(a::Number, f::CurlCurlGreen) = CurlCurlGreen(f.wavenumber, a*f.source, f.position)
