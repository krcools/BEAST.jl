struct PlaneWaveHH3DTD{T,P,F} <: TDFunctional{T}
    direction::P
    speed_of_light::T
    signature::F
    amplitude::T
end


function planewave(direction, speedoflight::Number, signature, amplitude=one(speedoflight)::Number)
    PlaneWaveHH3DTD(direction, speedoflight, signature, amplitude)
end

scalartype(p::PlaneWaveHH3DTD) = scalartype(p.amplitude)

*(a, f::PlaneWaveHH3DTD) = PlaneWaveHH3DTD(f.direction, f.speed_of_light, f.signature, a*f.amplitude)


function(f::PlaneWaveHH3DTD)(r,t)
    r = cartesian(r)
    t = cartesian(t)[1]

    k = f.direction
    u = sum(k[i]*r[i] for i in 1:length(k))
    #u = dot(f.direction, r)
    a = f.amplitude
    c = f.speed_of_light
    h = f.signature
    a * h(c*t - u)
end


function gradient(f::PlaneWaveHH3DTD)
    @assert f.amplitude ≈ 1
    PlaneWaveMWTD(
        f.direction,
        -f.direction,
        f.speed_of_light,
        derive(f.signature)
    )
end


struct DotTraceHH{T,F} <: TDFunctional{T}
    field::F
end

scalartype(d::DotTraceHH{T}) where {T} = T

dot(::NormalVector, field::TDFunctional) = DotTraceHH{scalartype(field), typeof(field)}(field)


function (f::DotTraceHH)(p,t)
    n = normal(p)
    x = cartesian(p)
    return dot(n, f.field(x,t))
end
