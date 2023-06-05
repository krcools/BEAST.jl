mutable struct PlaneWaveVIE{T,P} <: Functional
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
