struct NDLCDRefSpace{T} <: RefSpace{T} end

function valuetype(ref::NDLCDRefSpace{T}, charttype::Type) where {T}
    SVector{universedimension(charttype),T}
end

function (ϕ::NDLCDRefSpace)(ndlc)

    u,v,w = parametric(ndlc)
    j = jacobian(ndlc)

    tu = tangents(ndlc,1)
    tv = tangents(ndlc,2)
    tw = tangents(ndlc,3)

    # B = [tu tv tw]

    # return SVector((
    #     (value=(2*B*[(u-1),v,w]/j),divergence=(6/j)),
    #     (value=(2*B*[u,(v-1),w]/j),divergence=(6/j)),
    #     (value=(2*B*[u,v,(w-1)]/j),divergence=(6/j)),
    #     (value=(2*B*[u,v,w]/j)    ,divergence=(6/j))
    # ))

    return SVector((
        (value=(2*((u-1)*tu + v*tv + w*tw)/j), divergence=(6/j)),
        (value=(2*(u*tu + (v-1)*tv + w*tw)/j), divergence=(6/j)),
        (value=(2*(u*tu + v*tv + (w-1)*tw)/j),divergence=(6/j)),
        (value=(2*(u*tu + v*tv + w*tw)/j)    ,divergence=(6/j))
    ))
end

numfunctions(x::NDLCDRefSpace, dom::CompScienceMeshes.ReferenceSimplex{3}) = 4

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
    # el: the supporting chart of the input
    # q: the local index of the face of el on which we compute the trace
    # fc: the chart on which we compute the trace
    t = zeros(scalartype(x),3,4)
    fa = SArray{Tuple{2},Int,1,2}
    te = SArray{Tuple{3},Int,1,3}
    j = 0
    y = BEAST.NDRefSpace{scalartype(x)}()
    n = normal(fc)
    for fac in faces(el)
        j += 1
        j == q && continue
        # fac != fc || continue
        #vind locale index (i) in fc van de edge die in zowel fc als fac ligt
        A = setdiff(fc.vertices,setdiff(fc.vertices,fac.vertices))
        # @show length(A)
        @assert length(A) == 2
        edg = fa([findfirst(isequal(A[1]), fc.vertices),findfirst(isequal(A[2]), fc.vertices)])
        i = abs(CompScienceMeshes.relorientation(edg,te([1,2,3])))
        #i = findfirst(isequal(setdiff(fcv,facv)[1]), fcv)
        edg_chart = simplex(fc.vertices[edg[1]], fc.vertices[edg[2]])
        edg_centr = cartesian(center(edg_chart))
        p = neighborhood(fc, carttobary(fc, edg_centr))
        r = neighborhood(el, carttobary(el, edg_centr))
        yp = y(p)[i].value
        xr = x(r)[j].value
        tgt = edg_chart.vertices[1] - edg_chart.vertices[2]
        ypt = dot(yp, tgt)
        xrt = dot(n × xr, tgt)
        # t[i,j] = volume(el)/(volume(fac)*norm(A[1]-A[2]))
        t[i,j] = xrt / ypt
    end
    return t
end

divergence(ref::NDLCDRefSpace, sh, el) = [Shape(sh.cellid, 1, sh.coeff/volume(el))]


function restrict(ϕ::NDLCDRefSpace{T}, dom1, dom2) where {T}
    # dom2 is the smaller of the domains

    K = numfunctions(ϕ, domain(dom1))
    D = dimension(dom1)

    @assert K == 4
    @assert D == 3
    @assert D == dimension(dom2)

    Q = zeros(T,K,K)
    for (i,face) in enumerate(faces(dom2))

        p = center(face)
        c = cartesian(p)
        A = volume(face)

        m = normal(p)
        # u = carttobary(dom1,cartesian(p))

        u = carttobary(dom1, c)
        x = neighborhood(dom1, u)

        y = ϕ(x)

        for j in 1:K
            Q[j,i] = dot(y[j].value, m) * A
        end
    end

    return Q
end
