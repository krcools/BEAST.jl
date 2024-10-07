struct GWPCurlRefSpace{T,Degree} <: RefSpace{T} end

function numfunctions(x::GWPCurlRefSpace{<:Any,D},
        dom::CompScienceMeshes.ReferenceSimplex{2}) where{D} (D+1)*(D+3) end

function (ϕ::GWPCurlRefSpace{T,Degree})(p) where {T,Degree}
    dom = domain(chart(p))
    u = parametric(p)
    vals = ϕ(dom, u)
    pushforwardcurl(vals, p)
end

function (ϕ::GWPCurlRefSpace{T,Deg})(dom::CompScienceMeshes.ReferenceSimplex{Dim},
    u) where {T,Deg,Dim}

    d = Deg
    u, v = u
    w = 1-u-v

    s = range(zero(T), one(T), length=d+3)

    nd1 = point(T, -v, u-1)
    nd2 = point(T, -v+1, u)
    nd3 = point(T, -v, u)

    P = SVector{2,T}
    vals = P[]
    crls = T[]

    i = 0
    for j in 1:d+1
        k = (d+2)-i-j
            Rᵢ = _sylpoly(s, i+1, u)
            Rⱼ = _sylpoly_shift(s, j+1, v)
            Rₖ = _sylpoly_shift(s, k+1, w)
            push!(vals, Rᵢ*Rⱼ*Rₖ*nd1)

            dRᵢ = _sylpoly_diff(s, i+1, u)
            dRⱼ = _sylpoly_shift_diff(s, j+1, v)
            dRₖ = _sylpoly_shift_diff(s, k+1, w)
            du = dRᵢ*Rⱼ*Rₖ - Rᵢ*Rⱼ*dRₖ
            dv = Rᵢ*dRⱼ*Rₖ - Rᵢ*Rⱼ*dRₖ
            curl = (du*nd1[2] - dv*nd1[1]) + 2*Rᵢ*Rⱼ*Rₖ
            push!(crls, curl)
    end

    for i in 1:d+1
        j = 0
        k = (d+2)-i-j
        Rᵢ = _sylpoly_shift(s, i+1, u)
        Rⱼ = _sylpoly(s, j+1, v)
        Rₖ = _sylpoly_shift(s, k+1, w)
        push!(vals, Rᵢ*Rⱼ*Rₖ*nd2)

        dRᵢ = _sylpoly_shift_diff(s, i+1, u)
        dRⱼ = _sylpoly_diff(s, j+1, v)
        dRₖ = _sylpoly_shift_diff(s, k+1, w)
        du = dRᵢ*Rⱼ*Rₖ - Rᵢ*Rⱼ*dRₖ
        dv = Rᵢ*dRⱼ*Rₖ - Rᵢ*Rⱼ*dRₖ
        curl = (du*nd2[2] - dv*nd2[1]) + 2*Rᵢ*Rⱼ*Rₖ
        push!(crls, curl)
    end

    for i in 1:d+1
        j = (d+2)-i
        k = 0
        Rᵢ = _sylpoly_shift(s, i+1, u)
        Rⱼ = _sylpoly_shift(s, j+1, v)
        Rₖ = _sylpoly(s, k+1, w)
        push!(vals, Rᵢ*Rⱼ*Rₖ*nd3)

        dRᵢ = _sylpoly_shift_diff(s, i+1, u)
        dRⱼ = _sylpoly_shift_diff(s, j+1, v)
        dRₖ = _sylpoly_diff(s, k+1, w)

        du = dRᵢ*Rⱼ*Rₖ - Rᵢ*Rⱼ*dRₖ
        dv = Rᵢ*dRⱼ*Rₖ - Rᵢ*Rⱼ*dRₖ
        curl = (du*nd3[2] - dv*nd3[1]) + 2*Rᵢ*Rⱼ*Rₖ
        
        push!(crls, curl)
    end

    for i in 1:d+1
        for j in 1:d+1
            k = (d+2)-i-j
            k <= 0 && continue
            Rsᵢ = _sylpoly_shift(s, i+1, u)
            Rsⱼ = _sylpoly_shift(s, j+1, v)
            Rsₖ = _sylpoly_shift(s, k+1, w)
            Rᵢ = _sylpoly(s, i+1, u)
            Rⱼ = _sylpoly(s, j+1, v)
            Rₖ = _sylpoly(s, k+1, w)
            S1 = Rᵢ*Rsⱼ*Rsₖ*nd1
            S2 = Rsᵢ*Rⱼ*Rsₖ*nd2
            S3 = Rsᵢ*Rsⱼ*Rₖ*nd3
            push!(vals, S2-S3)
            push!(vals, S3-S1)

            dRsᵢ = _sylpoly_shift_diff(s, i+1, u)
            dRsⱼ = _sylpoly_shift_diff(s, j+1, v)
            dRsₖ = _sylpoly_shift_diff(s, k+1, w)
            dRᵢ = _sylpoly_diff(s, i+1, u)
            dRⱼ = _sylpoly_diff(s, j+1, v)
            dRₖ = _sylpoly_diff(s, k+1, w)

            du = dRᵢ*Rsⱼ*Rsₖ - Rᵢ*Rsⱼ*dRsₖ
            dv = Rᵢ*dRsⱼ*Rsₖ - Rᵢ*Rsⱼ*dRsₖ
            curlS1 = du*nd1[2] - dv*nd1[1] + 2*Rᵢ*Rsⱼ*Rsₖ

            du = dRsᵢ*Rⱼ*Rsₖ - Rsᵢ*Rⱼ*dRsₖ
            dv = Rsᵢ*dRⱼ*Rsₖ - Rsᵢ*Rⱼ*dRsₖ
            curlS2 = du*nd2[2] - dv*nd2[1] + 2*Rsᵢ*Rⱼ*Rsₖ

            du = dRsᵢ*Rsⱼ*Rₖ - Rsᵢ*Rsⱼ*dRₖ
            dv = Rsᵢ*dRsⱼ*Rₖ - Rsᵢ*Rsⱼ*dRₖ
            curlS3 = du*nd3[2] - dv*nd3[1] + 2*Rsᵢ*Rsⱼ*Rₖ

            push!(crls, curlS2 - curlS3)
            push!(crls, curlS3 - curlS1)
    end end

    NF = length(vals)
    SVector{NF}([(value=f, curl=c) for (f,c) in zip(vals, crls)])
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
        @show cartesian(p_edge)
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
    if d > 1
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
