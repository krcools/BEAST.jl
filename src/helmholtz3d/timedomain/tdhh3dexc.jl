struct PlaneWaveHH3DTD{T,P,F} <: TDFunctional{T}
    direction::P
    speed_of_light::T
    signature::F
    amplitude::T
end


#speedoflight(f::PlaneWaveHH3DTD) = f.speed_of_light


function planewave(direction, speedoflight::Number, signature, amplitude=one(speedoflight)::Number)
    PlaneWaveHH3DTD(direction, speedoflight, signature, amplitude)
end


*(a, f::PlaneWaveHH3DTD) = PlaneWaveHH3DTD(f.direction, f.speed_of_light, f.signature, a*f.amplitude)


function(f::PlaneWaveHH3DTD)(r,t)
    k = f.direction
    u = sum(k[i]*r[i] for i in 1:length(k))
    #u = dot(f.direction, r)
    a = f.amplitude
    c = f.speed_of_light
    h = f.signature
    a * h(c*t - u)
end
