mutable struct PlaneWaveVIE{T,P} <: Functional{T}
    direction::P
    polarisation::P
    wavenumber::T
    amplitude::T
  end

  scalartype(x::PlaneWaveVIE{T,P}) where {T,P} = Complex{T}

  function PlaneWaveVIE(d,p,k,a = 1)
    T = promote_type(eltype(d), eltype(p), typeof(k), typeof(a))
    P = similar_type(typeof(d), T)
    PlaneWaveVIE{T,P}(d,p,k,a)
end

"""
  planewavevie(;
      direction    = error("missing arguement `direction`"),
      polarization = error("missing arguement `polarization`"),
      wavenumber   = error("missing arguement `wavenumber`"),
      amplitude    = 1,
      ) 

For volume integral equations
"""
planewavevie(;
    direction    = error("missing arguement `direction`"),
    polarization = error("missing arguement `polarization`"),
    wavenumber   = error("missing arguement `wavenumber`"),
    amplitude    = 1,
    ) = PlaneWaveVIE(direction, polarization, wavenumber, amplitude)

function (e::PlaneWaveVIE)(p)
    k = e.wavenumber
    d = e.direction
    u = e.polarisation
    a = e.amplitude
    x = p
    a * exp(-im * k * dot(d, x)) * u
  end

  function (e::PlaneWaveVIE)(p::CompScienceMeshes.MeshPointNM)
    k = e.wavenumber
    d = e.direction
    u = e.polarisation
    a = e.amplitude
    x = cartesian(p)
    a * exp(-im * k * dot(d, x)) * u
  end

  function curl(field::PlaneWaveVIE)
    k = field.wavenumber
    d = field.direction
    u = field.polarisation
    a = field.amplitude
    v = d × u
    b = -a * im * k
    PlaneWaveVIE(d, v, k, b)
  end

*(a::Number, e::PlaneWaveVIE) = PlaneWaveVIE(e.direction, e.polarisation, e.wavenumber, a*e.amplitude)

integrand(::PlaneWaveVIE, test_vals, field_val) = test_vals[1] ⋅ field_val


# Excitation for Lippmann Schwinger Volume Integral Equation
mutable struct LinearPotentialVIE{T,P} <: Functional{T}
  direction::P
  amplitude::T
end

scalartype(x::LinearPotentialVIE{T,P}) where {T,P} = T

function LinearPotentialVIE_(d,a = 1)
  T = promote_type(eltype(d), typeof(a))
  P = similar_type(typeof(d), T) #SVector{3,T}
  return LinearPotentialVIE{T,P}(d,a)
end

"""
  linearpotentialvie(;
    direction    = error("missing argument `direction`"),
    amplitude    = 1,
  )

Linear potential for volume integral equations.
"""
linearpotentialvie(;
  direction    = error("missing argument `direction`"),
  amplitude    = 1,
) = LinearPotentialVIE_(direction, amplitude)

function (e::LinearPotentialVIE)(p)
  d = e.direction
  a = e.amplitude
  x = cartesian(p)
  return a * dot(d, x)
end

*(a::Number, e::LinearPotentialVIE) = LinearPotentialVIE_(e.direction, a*e.amplitude)

integrand(::LinearPotentialVIE, test_vals, field_val) = test_vals[1] ⋅ field_val