struct GWPCurlRefSpace{T,Degree} <: RefSpace{T} end

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
    # rwg1 = point(T, u-1, v)
    # rwg2 = point(T, u, v-1)
    # rwg3 = point(T, u, v)
    nd1 = point(T, -v, u-1)
    nd2 = point(T, -v+1, u)
    nd3 = point(T, -v, u)

    P = SVector{2,T}
    vals = P[]
    dffu = T[]
    dffv = T[]

    i = 0
    for j in 1:d+1
        k = (d+2)-i-j
            Rᵢ = _sylpoly(s, i+1, u)
            Rⱼ = _sylpoly_shift(s, j+1, v)
            Rₖ = _sylpoly_shift(s, k+1, w)
            push!(vals, Rᵢ*Rⱼ*Rₖ*nd1)

            # du = _sylpoly_shift_diff(s, i+1, u)*Rⱼ*Rₖ
    end

    for i in 1:d+1
        j = 0
        k = (d+2)-i-j
        Rᵢ = _sylpoly_shift(s, i+1, u)
        Rⱼ = _sylpoly(s, j+1, v)
        Rₖ = _sylpoly_shift(s, k+1, w)
        push!(vals, Rᵢ*Rⱼ*Rₖ*nd2)
    end

    for i in 1:d+1
        j = (d+2)-i
        k = 0
        Rᵢ = _sylpoly_shift(s, i+1, u)
        Rⱼ = _sylpoly_shift(s, j+1, v)
        Rₖ = _sylpoly(s, k+1, w)
        push!(vals, Rᵢ*Rⱼ*Rₖ*nd3)
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
            N1 = (d+2)/(d+2-i)
            N2 = (d+2)/(d+2-j)
            N3 = (d+2)/(d+2-k)
            push!(vals, S2-S3)
            push!(vals, S3-S1)
    end end

    NF = length(vals)
    SVector{NF}([(value=f, curl=zero(T)) for f in vals])
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