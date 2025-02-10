mutable struct PlaneWaveMWTD{T,F,P} <: TDFunctional{T}
  direction::P
  polarisation::P
  speedoflight::T
  amplitude::F
end

"""
    planewave(polarisation,direction,amplitude,speedoflight)

Time-domain plane wave.
"""
function planewave(polarisation,direction,amplitude,speedoflight)
    PlaneWaveMWTD(direction,polarisation,speedoflight,amplitude)
end

planewave(;signature, polarization, direction, speedoflight) =
    PlaneWaveMWTD(direction, polarization, speedoflight, signature)

scalartype(p::PlaneWaveMWTD) = eltype(p.polarisation)

*(a, pw::PlaneWaveMWTD) = PlaneWaveMWTD(
    pw.direction,
    a * pw.polarisation,
    pw.speedoflight,
    pw.amplitude
)

cross(k, pw::PlaneWaveMWTD) = PlaneWaveMWTD(
    pw.direction,
    k × pw.polarisation,
    pw.speedoflight,
    pw.amplitude
)

function (f::PlaneWaveMWTD)(r,t)
    t = cartesian(t)[1]
    r = cartesian(r)
    dr = zero(typeof(t))
    for i in 1 : 3
        dr += r[i]*f.direction[i]
    end
    f.polarisation * f.amplitude(t - dr/f.speedoflight)
end

function integrate(f::BEAST.PlaneWaveMWTD)
    planewave(
        signature = integrate(f.amplitude),
        direction = f.direction,
        polarization = f.polarisation,
        speedoflight = f.speedoflight)
end

function differentiate(f::BEAST.PlaneWaveMWTD)
    planewave(
        signature = derive(f.amplitude),
        direction = f.direction,
        polarization = f.polarisation,
        speedoflight = f.speedoflight)
end
