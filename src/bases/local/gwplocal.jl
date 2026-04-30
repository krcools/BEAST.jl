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

@generated function (::GWPCurlRefSpace{T,Degree})(dom::CompScienceMeshes.ReferenceSimplex{Dim},
    bary, ::Val{NF}) where {T,Degree,Dim,NF}
     
    u = :(bary[1])
    v = :(bary[2])
    w = :(one(T) - $u - $v)

    #nodes for sylvester polynomials
    s = :()
    for i in 0:Degree+2
        s = :($s..., $i/($Degree+2))
    end

    nd1 = :(SVector{2}( -$v, $u-one(T) ))
    nd2 = :(SVector{2}( -$v+one(T), $u ))
    nd3 = :(SVector{2}( -$v, $u ))
   
    nts = :()

    i = 0
    for j in 1:Degree+1
        k = (Degree+2)-i-j
        Rᵢ = gen_sylpoly(s, i+1, u, T)
        Rⱼ = gen_sylpoly_shift(s, j+1, v, T)
        Rₖ = gen_sylpoly_shift(s, k+1, w, T)
        
        dRᵢ = gen_sylpoly_diff(s, i+1, u, T)
        dRⱼ = gen_sylpoly_shift_diff(s, j+1, v, T)
        dRₖ = gen_sylpoly_shift_diff(s, k+1, w, T)
        du = :( $dRᵢ*$Rⱼ*$Rₖ - $Rᵢ*$Rⱼ*$dRₖ )
        dv = :( $Rᵢ*$dRⱼ*$Rₖ - $Rᵢ*$Rⱼ*$dRₖ )
        curl = :( ($du*$nd1[2] - $dv*$nd1[1]) + 2*$Rᵢ*$Rⱼ*$Rₖ )

        nts = :($nts..., (value=$Rᵢ*$Rⱼ*$Rₖ*$nd1, curl=$curl))
    end

    for i in 1:Degree+1
        j = 0
        k = (Degree+2)-i-j
        Rᵢ = gen_sylpoly_shift(s, i+1, u, T)
        Rⱼ = gen_sylpoly(s, j+1, v, T)
        Rₖ = gen_sylpoly_shift(s, k+1, w, T)
        
        dRᵢ = gen_sylpoly_shift_diff(s, i+1, u, T)
        dRⱼ = gen_sylpoly_diff(s, j+1, v, T)
        dRₖ = gen_sylpoly_shift_diff(s, k+1, w, T)
        du = :( $dRᵢ*$Rⱼ*$Rₖ - $Rᵢ*$Rⱼ*$dRₖ )
        dv = :( $Rᵢ*$dRⱼ*$Rₖ - $Rᵢ*$Rⱼ*$dRₖ )
        curl = :( ($du*$nd2[2] - $dv*$nd2[1]) + 2*$Rᵢ*$Rⱼ*$Rₖ )

        nts = :($nts..., (value=$Rᵢ*$Rⱼ*$Rₖ*$nd2, curl=$curl))
    end

    for i in 1:Degree+1
        j = (Degree+2)-i
        k = 0
        Rᵢ = gen_sylpoly_shift(s, i+1, u, T)
        Rⱼ = gen_sylpoly_shift(s, j+1, v, T)
        Rₖ = gen_sylpoly(s, k+1, w, T)
        
        dRᵢ = gen_sylpoly_shift_diff(s, i+1, u, T)
        dRⱼ = gen_sylpoly_shift_diff(s, j+1, v, T)
        dRₖ = gen_sylpoly_diff(s, k+1, w, T)
        
        du = :( $dRᵢ*$Rⱼ*$Rₖ - $Rᵢ*$Rⱼ*$dRₖ )
        dv = :( $Rᵢ*$dRⱼ*$Rₖ - $Rᵢ*$Rⱼ*$dRₖ )
        curl = :( ($du*$nd3[2] - $dv*$nd3[1]) + 2*$Rᵢ*$Rⱼ*$Rₖ )

        nts = :($nts..., (value=$Rᵢ*$Rⱼ*$Rₖ*$nd3, curl=$curl))
    end

    for i in 1:Degree+1
        for j in 1:Degree+1
            k = (Degree+2)-i-j
            k <= 0 && continue
            Rsᵢ = gen_sylpoly_shift(s, i+1, u, T)
            Rsⱼ = gen_sylpoly_shift(s, j+1, v, T)
            Rsₖ = gen_sylpoly_shift(s, k+1, w, T)
            Rᵢ = gen_sylpoly(s, i+1, u, T)
            Rⱼ = gen_sylpoly(s, j+1, v, T)
            Rₖ = gen_sylpoly(s, k+1, w, T)
            S1 = :( $Rᵢ*$Rsⱼ*$Rsₖ*$nd1 )
            S2 = :( $Rsᵢ*$Rⱼ*$Rsₖ*$nd2 )
            S3 = :( $Rsᵢ*$Rsⱼ*$Rₖ*$nd3 )

            dRsᵢ = gen_sylpoly_shift_diff(s, i+1, u, T)
            dRsⱼ = gen_sylpoly_shift_diff(s, j+1, v, T)
            dRsₖ = gen_sylpoly_shift_diff(s, k+1, w, T)
            dRᵢ = gen_sylpoly_diff(s, i+1, u, T)
            dRⱼ = gen_sylpoly_diff(s, j+1, v, T)
            dRₖ = gen_sylpoly_diff(s, k+1, w, T)

            du = :( $dRᵢ*$Rsⱼ*$Rsₖ - $Rᵢ*$Rsⱼ*$dRsₖ )
            dv = :( $Rᵢ*$dRsⱼ*$Rsₖ - $Rᵢ*$Rsⱼ*$dRsₖ )
            curlS1 = :( $du*$nd1[2] - $dv*$nd1[1] + 2*$Rᵢ*$Rsⱼ*$Rsₖ )

            du = :( $dRsᵢ*$Rⱼ*$Rsₖ - $Rsᵢ*$Rⱼ*$dRsₖ )
            dv = :( $Rsᵢ*$dRⱼ*$Rsₖ - $Rsᵢ*$Rⱼ*$dRsₖ )
            curlS2 = :( $du*$nd2[2] - $dv*$nd2[1] + 2*$Rsᵢ*$Rⱼ*$Rsₖ )

            du = :( $dRsᵢ*$Rsⱼ*$Rₖ - $Rsᵢ*$Rsⱼ*$dRₖ )
            dv = :( $Rsᵢ*$dRsⱼ*$Rₖ - $Rsᵢ*$Rsⱼ*$dRₖ )
            curlS3 = :( $du*$nd3[2] - $dv*$nd3[1] + 2*$Rsᵢ*$Rsⱼ*$Rₖ )

            nts = :($nts..., (value=$S2-$S3, curl=$curlS2 - $curlS3))
            
            nts = :($nts..., (value=$S3-$S1, curl=$curlS3 - $curlS1))
    end end

    ex = :(SVector{NF}(()))

    for i in 1:NF
        push!(ex.args[2].args, :($nts[$i]))
    end

    return ex
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
