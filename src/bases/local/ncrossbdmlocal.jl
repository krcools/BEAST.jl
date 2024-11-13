struct NCrossBDMRefSpace{T} <: RefSpace{T} end

function (f::NCrossBDMRefSpace{T})(p) where T

    u,v = parametric(p)
    n = normal(p)
    tu = tangents(p,1)
    tv = tangents(p,2)

    j = jacobian(p)
    d = 1/j

    return @SVector[
        (value= n × (-v*tu+v*tv)/j, curl=d),
        (value= n × ((u+v-1)*tu)/j,   curl=d),
        (value= n × ((u+v-1)*tv) /j,   curl=d),
        (value= n × (u*tu-u*tv)/j,  curl=d),
        (value= n × (u*tu)/j,         curl=d),
        (value= n × (v*tv)/j,         curl=d),]
end

numfunctions(x::NCrossBDMRefSpace, dom::CompScienceMeshes.ReferenceSimplex{2}) = 6
