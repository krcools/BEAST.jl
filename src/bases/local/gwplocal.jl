struct GWPCurlRefSpace{T,Degree} <: RefSpace{T} end

function numfunctions(x::GWPCurlRefSpace{<:Any,D},
        dom::CompScienceMeshes.ReferenceSimplex{2}) where {D}
        (D+1)*(D+3)
end
function dimtype(x::GWPCurlRefSpace{<:Any,D},
    dom::CompScienceMeshes.ReferenceSimplex{2}) where {D}
    Val((D+1)*(D+3))
end

function (ϕ::GWPCurlRefSpace{T,Degree})(p) where {T,Degree}
    dom = domain(chart(p))
    u = parametric(p)
    vals = ϕ(dom, u)
    pushforwardcurl(vals, p)
end

function (ϕ::GWPCurlRefSpace{T,Deg})(dom::CompScienceMeshes.ReferenceSimplex{Dim},
    u) where {T,Deg,Dim}

    ϕ(dom, u, dimtype(ϕ,dom))
end


macro gwp_shapefunction_barycentric_index(degree)
    k = degree + 2
    LUT = SVector{(degree+1)*(degree+3)}(vcat(
        [(0, j, k-j) for j in 1:degree+1],
        [(i, 0, k-i) for i in 1:degree+1],
        [(i, k-i, 0) for i in 1:degree+1],
        [repeat([(i, j, k-i-j)], 2) for i in 1:degree for j in 1:degree if k-i-j ≥ 1]...
    ))
    esc(LUT)
end

localidx(i,::Val{0}) =  @gwp_shapefunction_barycentric_index(0)[i]
localidx(i,::Val{1}) =  @gwp_shapefunction_barycentric_index(1)[i]
localidx(i,::Val{2}) =  @gwp_shapefunction_barycentric_index(2)[i]
localidx(i,::Val{3}) =  @gwp_shapefunction_barycentric_index(3)[i]
localidx(i,::Val{4}) =  @gwp_shapefunction_barycentric_index(4)[i]
localidx(i,::Val{5}) =  @gwp_shapefunction_barycentric_index(5)[i]
localidx(i,::Val{6}) =  @gwp_shapefunction_barycentric_index(6)[i]


function shapefunction(idx::Int, bary, ::Val{Degree}) where {Degree}

    T = eltype(bary)

    u, v = bary[1], bary[2]
    w = one(T)-v-u

    s = SVector{Degree+3,T}([i/(Degree+2) for i in 0:Degree+2])

    nd1 = SVector(-v, u-one(T))
    nd2 = SVector(-v+one(T), u)
    nd3 = SVector(-v, u)
    
    i,j,k = localidx(idx, Val(Degree))
    
    Rsᵢ = BEAST._sylpoly_shift(s, i+1, u)
    Rsⱼ = BEAST._sylpoly_shift(s, j+1, v)
    Rsₖ = BEAST._sylpoly_shift(s, k+1, w)
    
    Rᵢ = BEAST._sylpoly(s, i+1, u)
    Rⱼ = BEAST._sylpoly(s, j+1, v)
    Rₖ = BEAST._sylpoly(s, k+1, w)
    
    S1 = Rᵢ*Rsⱼ*Rsₖ*nd1
    S2 = Rsᵢ*Rⱼ*Rsₖ*nd2
    S3 = Rsᵢ*Rsⱼ*Rₖ*nd3
    
    dRsᵢ = BEAST._sylpoly_shift_diff(s, i+1, u)
    dRsⱼ = BEAST._sylpoly_shift_diff(s, j+1, v)
    dRsₖ = BEAST._sylpoly_shift_diff(s, k+1, w)
    dRᵢ = BEAST._sylpoly_diff(s, i+1, u)
    dRⱼ = BEAST._sylpoly_diff(s, j+1, v)
    dRₖ = BEAST._sylpoly_diff(s, k+1, w)
    
    du = dRᵢ*Rsⱼ*Rsₖ - Rᵢ*Rsⱼ*dRsₖ
    dv = Rᵢ*dRsⱼ*Rsₖ - Rᵢ*Rsⱼ*dRsₖ
    curlS1 = du*nd1[2] - dv*nd1[1] + 2*Rᵢ*Rsⱼ*Rsₖ
    
    du = dRsᵢ*Rⱼ*Rsₖ - Rsᵢ*Rⱼ*dRsₖ
    dv = Rsᵢ*dRⱼ*Rsₖ - Rsᵢ*Rⱼ*dRsₖ
    curlS2 = du*nd2[2] - dv*nd2[1] + 2*Rsᵢ*Rⱼ*Rsₖ
    
    du = dRsᵢ*Rsⱼ*Rₖ - Rsᵢ*Rsⱼ*dRₖ
    dv = Rsᵢ*dRsⱼ*Rₖ - Rsᵢ*Rsⱼ*dRₖ
    curlS3 = du*nd3[2] - dv*nd3[1] + 2*Rsᵢ*Rsⱼ*Rₖ

    if idx <= (Degree+1) 
        val = S1
        curl = curlS1
    elseif idx <= 2*(Degree+1)
        val = S2
        curl = curlS2
    elseif idx <= 3*(Degree+1)
        val = S3
        curl = curlS3
    else
        if (idx-3*(Degree+1)) % 2 == 1
            val = S2 - S3
            curl = curlS2 - curlS3
        else 
            val = S3 - S1
            curl = curlS3 - curlS1
        end
    end

    return (value=val, curl=curl)
end

function (::GWPCurlRefSpace{T,Degree})(dom::CompScienceMeshes.ReferenceSimplex{Dim},
    bary, ::Val{NF}) where {T,Degree,Dim,NF}
     
    return SVector{NF}(shapefunction(i, bary, Val(Degree)) for i in 1:NF)
end


function interpolate(fields, interpolant::GWPCurlRefSpace{T,Degree}, chart) where {T,Degree}

    d = Degree
    dim = (d+1)*(d+3)

    s = range(zero(T), one(T), length=d+3)

    edges = faces(chart)

    edge = edges[1]
    fields_edge = trace(edge, chart, fields)
    i = 0
    Q1 = stack(1:d+1) do j
        k = (d+2)-i-j
        u_edge = s[j+1]
        p_edge = neighborhood(edge, (u_edge,))
        # @show cartesian(p_edge)
        t_edge = -tangents(p_edge, 1)
        vals = fields_edge(u_edge)
        [dot(t_edge, val) for val in vals]
    end

    edge = edges[2]
    fields_edge = trace(edge, chart, fields)
    Q2 = stack(1:d+1) do i
        j = 0
        k = (d+2)-i-j
        u_edge = 1 - s[i+1]
        p_edge = neighborhood(edge, (u_edge,))
        t_edge = -tangents(p_edge, 1)
        vals = fields_edge(u_edge)
        [dot(t_edge, val) for val in vals]
    end

    edge = edges[3]
    fields_edge = trace(edge, chart, fields)
    Q3 = stack(1:d+1) do i
        j = (d+2)-i
        k = 0
        u_edge = 1-s[j+1]
        p_edge = neighborhood(edge, (u_edge,))
        t_edge = -tangents(p_edge, 1)
        vals = fields_edge(u_edge)
        [dot(t_edge, val) for val in vals]
    end

    Q = hcat(Q1,Q2,Q3)
    if d >= 1
        S = ((i,j,d+2-i-j) for i in 1:d+1 for j in 1:d+1 if d+2-i-j > 0)
        for (i,j,k) in S
            p_chart = neighborhood(chart, (s[i+1],s[j+1]))
            t_i = tangents(p_chart, 1)
            t_j = tangents(p_chart, 2)
            vals = fields(p_chart)
            q_i = [dot(t_i, val) for val in vals]
            q_j = [dot(t_j, val) for val in vals]
            Q = hcat(Q, q_i, q_j)
        end
    end
    return Q
end

function localindices(localspace::GWPCurlRefSpace{<:Any,Degree}, domain,
    dim::Type{Val{1}}, i) where {Degree}
    
    ne = Degree+1
    (i-1)*ne .+ (1:ne)
end

function localindices(localspace::GWPCurlRefSpace{<:Any,Degree}, domain,
    dim::Type{Val{2}}, i) where {Degree}
    
    ne = Degree+1
    nf = Degree * (Degree + 1)
    3*ne .+ (1:nf)
end
