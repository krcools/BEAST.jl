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

mutable struct ElectricDipole{T,P} <: Dipole
    location::P
    orientation::P
    wavenumber::T
    ε::T
    μ::T
end

function ElectricDipole(l,o,k,ε,μ)
    T = promote_type(eltype(l), eltype(o), typeof(k), typeof(ε), typeof(μ))
    P = similar_type(typeof(l), T)
    ElectricDipole{T,P}(l,o,k,ε,μ)
end

mutable struct curlElectricDipole{T,P} <: Dipole
    location::P
    orientation::P
    wavenumber::T
    ε::T
    μ::T
end

function curlElectricDipole(l,o,k,ε,μ)
    T = promote_type(eltype(l), eltype(o), typeof(k), typeof(ε), typeof(μ))
    P = similar_type(typeof(l), T)
    curlElectricDipole{T,P}(l,o,k,ε,μ)
end

mutable struct MagneticDipole{T,P} <: Dipole
    location::P
    orientation::P
    wavenumber::T
    ε::T
    μ::T
end

function MagneticDipole(l,o,k,ε,μ)
    T = promote_type(eltype(l), eltype(o), typeof(k), typeof(ε), typeof(μ))
    P = similar_type(typeof(l), T)
    MagneticDipole{T,P}(l,o,k,ε,μ)
end

mutable struct curlMagneticDipole{T,P} <: Dipole
    location::P
    orientation::P
    wavenumber::T
    ε::T
    μ::T
end

function curlMagneticDipole(l,o,k,ε,μ)
    T = promote_type(eltype(l), eltype(o), typeof(k), typeof(ε), typeof(μ))
    P = similar_type(typeof(l), T)
    curlMagneticDipole{T,P}(l,o,k,ε,μ)
end

"""
    electricdipole(;location, orientation, wavenumber, permittivity, permeability)

Create an electric dipole solution to Maxwell's equations representing the electric
field part. Implementation is based on (9.18) of Jackson's “Classical electrodynamics”,
with the notable difference that the ``\exp(ikr)`` is used.
"""
electricdipole(;
    location    = error("missing arguement `location`"),
    orientation = error("missing arguement `orientation`"),
    wavenumber   = error("missing arguement `wavenumber`"),
    permittivity = error("missing arguement `permittivity`"),
    permeability = error("missing arguement `permeability`")) =
    ElectricDipole(location, orientation, wavenumber, permittivity, permeability)

function (hd::ElectricDipole)(x; isfarfield=false)
    k = hd.wavenumber
    x_0 = hd.location
    p = hd.orientation
    r = norm(x-x_0)
    n =  (x - x_0)/r
    η = sqrt(hd.μ/hd.ε)
    if isfarfield
      # postfactor (4*π*im)/k to be consistent with BEAST far field computation
      # and, of course, omitted exp(-im*k*r)/r factor in (9.19)
      return η*k^2/(4*π*sqrt(hd.ε*hd.μ))*cross(cross(n,p),n)*(4*π*im)/k
    else
      return (1/(4*π*hd.ε))*exp(-im*k*r)*(k^2/r*cross(cross(n,p),n) + (1/r^3 + im*k/r^2)*(3*n*dot(n,p) - p))
    end
end

function (hd::curlElectricDipole)(x; isfarfield=false)
    k = hd.wavenumber
    x_0 = hd.location
    p = hd.orientation
    r = norm(x-x_0)
    n =  (x - x_0)/r
    c = 1/sqrt(hd.ε*hd.μ)
    if isfarfield
      # prefactor (-im*hd.μ*c*k), because this is the curl
      # postfactor (4*π*im)/k to be consistent with BEAST far field computation
      return (-im*hd.μ*c*k)*k^2/(4*π*sqrt(hd.ε*hd.μ))*cross(n,p)*(4*π*im)/k
    else
      #return (k^2/(4*π*sqrt(hd.ε*hd.μ)))*cross(n,p)*exp(-im*k*r)/r*(1 + 1/(im*k*r))
      return -im*hd.μ*c^2*(k^3)/(4*π)*cross(n,p)*exp(-im*k*r)/r*(1 + 1/(im*k*r))
    end
end

function curl(ehd::ElectricDipole)
    return curlElectricDipole(ehd.location,ehd.orientation,ehd.wavenumber,ehd.ε,ehd.μ)
end

"""
    magneticdipole(;location, orientation, wavenumber, permittivity, permeability)

Create a magnetic dipole solution to Maxwell's equations representing the electric
field part. Implementation is based on (9.36) of Jackson's “Classical electrodynamics”,
with the notable difference that the ``\exp(ikr)`` is used.
"""
magneticdipole(;
    location    = error("missing arguement `location`"),
    orientation = error("missing arguement `orientation`"),
    wavenumber   = error("missing arguement `wavenumber`"),
    permittivity = error("missing arguement `permittivity`"),
    permeability = error("missing arguement `permeability`")) =
    MagneticDipole(location, orientation, wavenumber, permittivity, permeability)

function (hd::MagneticDipole)(x; isfarfield=false)
    k = hd.wavenumber
    x_0 = hd.location
    m = hd.orientation
    r = norm(x-x_0)
    n =  (x - x_0)/r
    η = sqrt(hd.μ/hd.ε)
    if isfarfield
      # Jackson (9.36) without exp(-im*k*r)/r, r → ∞ and with (im*4*π)/k to match BEAST's farfield
      return -η/(4*π)*k^2*cross(n,m)*(im*4*π)/k 
    else
      # Jackson (9.36)
      return -η/(4*π)*k^2*cross(n,m)*exp(-im*k*r)/r*(1 + 1/(im*k*r))
    end
end

function (hd::curlMagneticDipole)(x; isfarfield=false)
    k = hd.wavenumber
    x_0 = hd.location
    m = hd.orientation
    r = norm(x-x_0)
    n =  (x - x_0)/r
    c = 1/sqrt(hd.ε*hd.μ)
    if isfarfield
       # Jackson (9.35) without exp(-im*k*r)/r, r → ∞ and with (im*4*π)/k to match BEAST's farfield
      return -im*hd.μ*c*k/(4π)*(k^2*cross(cross(n,m),n))*(im*4*π)/k 
    else
      # Jackson (9.35)
      return -im*hd.μ*c*k/(4π)*exp(-im*k*r)*(k^2*cross(cross(n,m),n)/r + (3*n*dot(n,m)-m)*(1/r^3 + im*k/r^2))
    end
end

function curl(ehd::MagneticDipole)
    return curlMagneticDipole(ehd.location,ehd.orientation,ehd.wavenumber,ehd.ε,ehd.μ)
end

*(a::Number, e::ElectricDipole) = ElectricDipole(e.location, a .* e.orientation, e.wavenumber,e.ε,e.μ)
*(a::Number, e::curlElectricDipole) = curlElectricDipole(e.location, a .* e.orientation, e.wavenumber,e.ε,e.μ)
*(a::Number, e::MagneticDipole) = MagneticDipole(e.location, a .* e.orientation, e.wavenumber,e.ε,e.μ)
*(a::Number, e::curlMagneticDipole) = curlMagneticDipole(e.location, a .* e.orientation, e.wavenumber,e.ε,e.μ)

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
