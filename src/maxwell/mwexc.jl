mutable struct PlaneWaveMW{T,P}
  direction::P
  polarisation::P
  gamma::T
  amplitude::T
end

function PlaneWaveMW(d,p,γ,a = 1)
    T = promote_type(eltype(d), eltype(p), typeof(γ), typeof(a))
    P = similar_type(typeof(d), T)
    PlaneWaveMW{T,P}(d,p,γ,a)
end

scalartype(x::PlaneWaveMW{T,P}) where {T,P} = promote_type(T, eltype(P)) 

"""
    planewavemw3d(;direction, polarization, wavenumber, gamma[, amplitude=1])

Create a plane wave solution to Maxwell's equations.
"""
function planewavemw3d(;
    direction    = error("missing arguement `direction`"),
    polarization = error("missing arguement `polarization`"),
    wavenumber   = nothing,
    gamma        = nothing,
    amplitude    = 1,
    )

    gamma, wavenumber = gamma_wavenumber_handler(gamma, wavenumber)
    isstatic(gamma) && (gamma = zero(eltype(direction)))

    return PlaneWaveMW(direction, polarization, gamma, amplitude)

end

function (e::PlaneWaveMW)(x)
  γ = e.gamma
  d = e.direction
  u = e.polarisation
  a = e.amplitude
  a * exp(-γ * dot(d, x)) * u
end

function curl(field::PlaneWaveMW)
  γ = field.gamma
  d = field.direction
  u = field.polarisation
  a = field.amplitude
  v = d × u
  b = -a * γ
  PlaneWaveMW(d, v, γ, b)
end

*(a::Number, e::PlaneWaveMW) = PlaneWaveMW(e.direction, e.polarisation, e.gamma, a*e.amplitude)

abstract type Dipole end

mutable struct DipoleMW{T,P} <: Dipole
    location::P
    orientation::P
    gamma::T
end

function DipoleMW(l,o,γ)
    T = promote_type(eltype(l), eltype(o), typeof(γ))
    P = similar_type(typeof(l), T)
    DipoleMW{T,P}(l,o,γ)
end

scalartype(x::DipoleMW{T,P}) where {T,P} = promote_type(T, eltype(P)) 

mutable struct curlDipoleMW{T,P} <: Dipole
    location::P
    orientation::P
    gamma::T
end

function curlDipoleMW(l,o,γ)
    T = promote_type(eltype(l), eltype(o), typeof(γ))
    P = similar_type(typeof(l), T)
    curlDipoleMW{T,P}(l,o,γ)
end

scalartype(x::curlDipoleMW{T,P}) where {T,P} = promote_type(T, eltype(P)) 

"""
    dipolemw3d(;location, orientation, wavenumber)

Create an electric dipole solution to Maxwell's equations representing the electric
field part. Implementation is based on (9.18) of Jackson's “Classical electrodynamics”,
with the notable difference that the ``\exp(ikr)`` is used.
"""
function dipolemw3d(;
    location    = error("missing arguement `location`"),
    orientation = error("missing arguement `orientation`"),
    wavenumber  = nothing,
    gamma       = nothing
    )

    gamma, wavenumber = gamma_wavenumber_handler(gamma, wavenumber)
    isstatic(gamma) && (gamma = zero(eltype(orientation)))

    return DipoleMW(location, orientation, gamma)

end

function (d::DipoleMW)(x; isfarfield=false)
    γ = d.gamma
    x_0 = d.location
    p = d.orientation
    if isfarfield
      # postfactor (4*π*im)/k = (-4*π)/γ to be consistent with BEAST far field computation
      # and, of course, adapted phase factor exp(im*k*dot(n,x_0)) with
      # respect to (9.19) of Jackson's Classical Electrodynamics
      r = norm(x)
      n = x/r
      return cross(cross(n,1/(4*π)*(-γ^2*cross(cross(n,p),n))*exp(γ*dot(n,x_0))*(-4*π)/γ),n)
    else
      r = norm(x-x_0)
      n = (x - x_0)/r
      return 1/(4*π)*exp(-γ*r)*(-γ^2/r*cross(cross(n,p),n) + 
              (1/r^3 + γ/r^2)*(3*n*dot(n,p) - p))
    end
end

function (d::curlDipoleMW)(x; isfarfield=false)
    γ = d.gamma
    x_0 = d.location
    p = d.orientation
    if isfarfield
      # postfactor (4*π*im)/k to be consistent with BEAST far field computation
      r = norm(x)
      n = x/r
      return (-γ)*(-γ^2)/(4*π)*cross(n,p)*(-4*π)/γ*exp(γ*dot(n,x_0))
    else
      r = norm(x-x_0)
      n =  (x - x_0)/r
      return γ^3/(4*π)*cross(n,p)*exp(-γ*r)/r*(1 + 1/(γ*r))
    end
end

function curl(d::DipoleMW)
    return curlDipoleMW(d.location, d.orientation, d.gamma)
end

*(a::Number, d::DipoleMW) = DipoleMW(d.location, a .* d.orientation, d.gamma)
*(a::Number, d::curlDipoleMW) = curlDipoleMW(d.location, a .* d.orientation, d.gamma)

mutable struct CrossTraceMW{F} <: Functional
  field::F
end

scalartype(x::CrossTraceMW) = scalartype(x.field)

mutable struct TangTraceMW{F} <: Functional
  field::F
end

scalartype(t::TangTraceMW) = scalartype(t.field)

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

struct NDotTrace{T,F} <: Functional
  field::F
end

NDotTrace(f::F) where {F} = NDotTrace{scalartype(f), F}(f)
NDotTrace{T}(f::F) where {T,F} = NDotTrace{T,F}(f)
scalartype(s::NDotTrace{T}) where {T} = T

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
