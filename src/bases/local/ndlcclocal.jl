struct NDLCCRefSpace{T} <: RefSpace{T} end

function valuetype(ref::NDLCCRefSpace{T}, charttype::Type) where {T}
    SVector{universedimension(charttype),T}
end

function (ϕ::NDLCCRefSpace)(ndlc)

    u,v,w = parametric(ndlc)
    j = jacobian(ndlc)

    tu = tangents(ndlc,1)
    tv = tangents(ndlc,2)
    tw = tangents(ndlc,3)

    #B = [tv-tu tw-tu -tu]
    B = [tu tv tw]
    BMT = transpose(inv(B))

    du = SVector{3}(BMT[:,1])
    dv = SVector{3}(BMT[:,2])
    dw = SVector{3}(BMT[:,3])

    #Choose 1-(0,0,0), 2-(1,0,0), 3-(0,1,0), 4-(0,0,1)
    #Now it is listed as:
    #(1,2) [1-v-w,u,u]  (direction: from 1 to 2)
    #(1,3) [v,1-u-w,v]
    #(1,4) [w,w,1-u,v]
    #(2,3) [-v,u,0]
    #(2,4) [-w,0,u]
    #(3,4) [0,-w,v]

    # return SVector((
    #     (value=(BMT*[-v,u,0])       ,curl=(B*[0,0,2])/j) ,
    #     (value=(BMT*[-w,0,u])       ,curl=(B*[0,-2,0])/j),
    #     (value=(BMT*[(w+v-1),-u,-u]),curl=(B*[0,2,-2])/j),
    #     (value=(BMT*[0,-w,v])       ,curl=(B*[2,0,0])/j) ,
    #     (value=(BMT*[-v,(u+w-1),-v]),curl=(B*[-2,0,2])/j),
    #     (value=(BMT*[-w,-w,(u+v-1)]),curl=(B*[2,-2,0])/j)
    # ))

    return SVector((
        (value=-v*du + u*dv, curl=(2*tw)/j),
        (value=-w*du + u*dw, curl=(-2*tv)/j),
        (value=(w+v-1)*du - u*dv - u*dw, curl=(2*tv - 2*tw)/j),
        (value=-w*dv + v*dw, curl=(2*tu)/j),
        (value=-v*du + (u+w-1)*dv -v*dw, curl=(-2*tu + 2*tw)/j),
        (value=-w*du -w*dv + (u+v-1)*dw, curl=(2*tu -2*tv)/j),
    ))

##    	return SVector((
##        (value=(BMT*[(w+v-1),-u,-u]),curl=(B*[0,2,-2])/j),
##        (value=(BMT*[-v,(u+w-1),-v])  ,curl=(B*[-2,0,2])/j),
##        (value=(BMT*[-w,-w,(u+v-1)])  ,curl=(B*[2,-2,0])/j),
##        (value=(BMT*[-v,u,0])       ,curl=(B*[0,0,2])/j) ,
##        (value=(BMT*[-w,0,u])       ,curl=(B*[0,-2,0])/j),
##        (value=(BMT*[0,-w,v])       ,curl=(B*[2,0,0])/j)
##        ))
end

numfunctions(x::NDLCCRefSpace, dom::CompScienceMeshes.ReferenceSimplex{3}) = 6

#check orientation
function curl(ref::NDLCCRefSpace, sh, el)
    a = [4,2,3,4,1,2]##[2,1,1,1,2,4]#[2,1,4,4,3,2]#[4,2,3,4,1,2]
    b = [3,4,2,1,3,1]##[3,2,4,3,4,3]#[4,3,3,1,2,1]#[3,4,2,1,3,1]
    sh1 = Shape(sh.cellid, b[sh.refid], -sh.coeff)
    sh2 = Shape(sh.cellid, a[sh.refid], sh.coeff)
    return [sh1,sh2]
end


function restrict(ϕ::NDLCCRefSpace{T}, dom1, dom2) where {T}
    # dom2 is the smaller of the domains

    K = numfunctions(ϕ, domain(dom1))
    D = dimension(dom1)

    @assert K == 6
    @assert D == 3
    @assert D == dimension(dom2)

    Q = zeros(T,K,K)
    for (i,edge) in enumerate(CompScienceMeshes.edges(dom2))

        p = center(edge)
        c = cartesian(p)
        A = volume(edge)

        t = tangents(p,1)
        t = normalize(t)
        u = carttobary(dom1, c)
        x = neighborhood(dom1, u)

        y = ϕ(x)

        for j in 1:K
            Q[j,i] = -dot(y[j].value, t) * A
        end
    end

    return Q
end


function ttrace(x::NDLCCRefSpace, el, q, fc)
    # el: the supporting chart of the input
    # q: the local index of the face of el on which we compute the trace
    # fc: the chart on which we compute the trace

    el_ctr = cartesian(center(el))
    fc_ctr = cartesian(center(fc))

    T = scalartype(x)
    y = BEAST.RTRefSpace{T}()
    t = zeros(scalartype(x),3,6)
    for (j,edgej) in enumerate(CompScienceMeshes.edges(el))
        for (i,edgei) in enumerate(CompScienceMeshes.edges(fc))
            nbdi = center(edgei)
            nbdj = neighborhood(el, carttobary(el, cartesian(nbdi)))
            @assert isapprox(cartesian(nbdi), cartesian(nbdj), atol=1e-4)
            xval = x(nbdj)[j].value
            # yval = y(nbdi)[i].value
            tgt = -tangents(nbdi,1)
            nrm = normal(fc)
            bnr = cross(normalize(tgt), nrm)
            @assert dot(bnr, cartesian(nbdi)-fc_ctr) > 0
            @assert norm(bnr) ≈ 1
            dot(nrm, fc_ctr - el_ctr) < 0 && (nrm = -nrm)
            t[i,j] = dot(bnr, nrm × xval) * volume(edgei)
        end
    end
    return t
    #
    # # fa = SArray{Tuple{2},Int,1,2}
    # # te = SArray{Tuple{3},Int,1,3}
    # j = 0
    # y = BEAST.NDRefSpace{scalartype(x)}()
    # n = normal(fc)
    # for fac in faces(el)
    #     j += 1
    #     j == q && continue
    #     # fac != fc || continue
    #     #vind locale index (i) in fc van de edge die in zowel fc als fac ligt
    #     A = setdiff(fc.vertices,setdiff(fc.vertices,fac.vertices))
    #     # @show length(A)
    #     @assert length(A) == 2
    #     edg = fa([findfirst(isequal(A[1]), fc.vertices),findfirst(isequal(A[2]), fc.vertices)])
    #     i = abs(CompScienceMeshes.relorientation(edg,te([1,2,3])))
    #     #i = findfirst(isequal(setdiff(fcv,facv)[1]), fcv)
    #     edg_chart = simplex(fc.vertices[edg[1]], fc.vertices[edg[2]])
    #     edg_centr = cartesian(center(edg_chart))
    #     p = neighborhood(fc, carttobary(fc, edg_centr))
    #     r = neighborhood(el, carttobary(el, edg_centr))
    #     yp = y(p)[i].value
    #     xr = x(r)[j].value
    #     tgt = edg_chart.vertices[1] - edg_chart.vertices[2]
    #     ypt = dot(yp, tgt)
    #     xrt = dot(n × xr, tgt)
    #     # t[i,j] = volume(el)/(volume(fac)*norm(A[1]-A[2]))
    #     t[i,j] = xrt / ypt
    # end
    # return t
end
