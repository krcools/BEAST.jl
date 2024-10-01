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
            Rᵢ = _sylpoly_shift(s, i+1, u)
            Rⱼ = _sylpoly(s, j+1, v)
            Rₖ = _sylpoly(s, k+1, w)
            push!(vals, Rᵢ*Rⱼ*Rₖ*nd1)

            # du = _sylpoly_shift_diff(s, i+1, u)*Rⱼ*Rₖ
    end

    for i in 1:d+1
        j = 0
        k = (d+2)-i-j
        Rᵢ = _sylpoly(s, i+1, u)
        Rⱼ = _sylpoly_shift(s, j+1, v)
        Rₖ = _sylpoly(s, k+1, w)
        push!(vals, Rᵢ*Rⱼ*Rₖ*nd2)
    end

    for i in 1:d+1
        j = (d+2)-i
        k = 0
        Rᵢ = _sylpoly(s, i+1, u)
        Rⱼ = _sylpoly(s, j+1, v)
        Rₖ = _sylpoly_shift(s, k+1, w)
        push!(vals, Rᵢ*Rⱼ*Rₖ*nd3)
    end

    for i in 1:d+1
        for j in 1:d+1
            k = (d+2)-i-j
            k <= 0 && continue
            Rsᵢ = _sylpoly_shift(s, i+1, u)
            Rsⱼ = _sylpoly_shift(s, j+1, v)
            Rᵢ = _sylpoly(s, i+1, u)
            Rⱼ = _sylpoly(s, j+1, v)
            Rₖ = _sylpoly(s, k+1, w)
            push!(vals, Rsᵢ*Rⱼ*Rₖ*nd1)
            push!(vals, Rᵢ*Rsⱼ*Rₖ*nd2)
    end end

    NF = length(vals)
    SVector{NF}([(value=f, curl=zero(T)) for f in vals])
end


