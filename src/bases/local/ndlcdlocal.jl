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
        (value=(2*B*[(u-1),v,w]/j),divergence=(6/j)),
        (value=(2*B*[u,(v-1),w]/j),divergence=(6/j)),
        (value=(2*B*[u,v,(w-1)]/j),divergence=(6/j)),
        (value=(2*B*[u,v,w]/j)    ,divergence=(6/j))
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
    fa = SArray{Tuple{2},Int,1,2}
    te = SArray{Tuple{3},Int,1,3}
    j = 0
    for fac in faces(el)
        j += 1
        fac != fc || continue
        #vind locale index (i) in fc van de edge die in zowel fc als fac ligt
        A = setdiff(fc.vertices,setdiff(fc.vertices,fac.vertices))
        edg = fa([findfirst(isequal(A[1]), fc.vertices),findfirst(isequal(A[2]), fc.vertices)])
        i = abs(CompScienceMeshes.relorientation(edg,te([1,2,3])))
        #i = findfirst(isequal(setdiff(fcv,facv)[1]), fcv)
        t[i,j] = volume(el)/(volume(fac)*norm(A[1]-A[2]))
    end
    return t
end
