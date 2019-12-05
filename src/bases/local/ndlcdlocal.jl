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

    return SVector((
        (value=(B*[(u-1),v,w]/j),divergence=(3/j)),
        (value=(B*[u,(v-1),w]/j),divergence=(3/j)),
        (value=(B*[u,v,(w-1)]/j),divergence=(3/j)),
        (value=(B*[u,v,w]/j)    ,divergence=(3/j))
    ))
end

function ntrace(x::NDLCDRefSpace, el, q, fc)
    t = zeros(scalartype(x),1,4)
    t[q] = 1 / volume(fc)
    return t
end

"""
Does not give the correct result for an imput basis with non-vanishing
input basis on the boundary
"""
function ttrace(x::NDLCDRefSpace, el, q, fc)
    t = zeros(scalartype(x),3,4)
    for fac in faces(el)
        fac != fc || continue
        #vind locale index (i) in fc van de edge die in zowel fc als fac ligt
        fcv = fc.vertices
        facv = fac.vertices
        i = findfirst(isequal(setdiff(fcv,facv)[1]), fcv)
        print(i)
        t[i,q] = 1
    end
    return t
end
