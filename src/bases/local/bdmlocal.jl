struct BDMRefSpace{T} <: RefSpace{T,6} end

function (f::BDMRefSpace)(p)

    u,v = parametric(p)

    tu = tangents(p,1)
    tv = tangents(p,2)

    j = jacobian(p)
    d = 1/j

    return @SVector[
        (value=(-v*tu+v*tv)/j, divergence=d),
        (value=(u+v-1)*tu/j,   divergence=d),
        (value=(u+v-1)*tv/j,   divergence=d),
        (value=(u*tu-u*tv)/j,  divergence=d),
        (value=u*tu/j,         divergence=d),
        (value=v*tv/j,         divergence=d),]
end
