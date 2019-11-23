struct NDLCDRefSpace{T} <: RefSpace{T,4} end

function valuetype(ref::NDLCDRefSpace{T}, charttype::Type) where {T}
    SVector{universedimension(charttype),T}
end

function (ϕ::NDLCDRefSpace)(ndlc)

    u,v,w = parametric(ndlc)
    j = jacobian(ndlc)

    tu = tangents(ndlc,1)
    tv = tangents(ndlc,2)
    tw = tangents(ndlc,3)

    B = [tu tv tw]

    #Choose 1-(0,0,0), 2-(1,0,0), 3-(0,1,0), 4-(0,0,1)
    #Now it is listed as:
    #(1,2,3) [u,v,w-1]
    #(1,2,4) [u,v-1,w]
    #(1,3,4) [u-1,v,w]
    #(2,3,4) [u,v,w]

    return SVector((
        (value=(B*[(u-1),v,w]/j),divergence=(3/j)),
        (value=(B*[u,(v-1),w]/j),divergence=(3/j)),
        (value=(B*[u,v,(w-1)]/j),divergence=(3/j)),
        (value=(B*[u,v,w]/j)    ,divergence=(3/j))
    ))
end
