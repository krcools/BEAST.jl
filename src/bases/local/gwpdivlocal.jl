struct GWPDivRefSpace{T,Degree} <: RefSpace{T} end

function numfunctions(x::GWPDivRefSpace{<:Any,D},
    dom::CompScienceMeshes.ReferenceSimplex{2}) where{D} (D+1)*(D+3) end

function (ϕ::GWPDivRefSpace{T,Degree})(p) where {T,Degree}
    ψ = GWPCurlRefSpace{T,Degree}()
    dom = domain(chart(p))
    u = parametric(p)
    vals = ψ(dom, u)
    vals = pushforwardcurl(vals, p)
    n = normal(p)
    map(vals) do v
        (value = -n × v.value, divergence=v.curl)
    end
end

function interpolate(fields, interpolant::GWPDivRefSpace{T,Degree}, chart) where {T,Degree}

    function nxfields(p)
        vals = fields(p)
        n = normal(p)
        return map(v -> n×v, vals)
    end

    return interpolate(nxfields, GWPCurlRefSpace{T,Degree}(), chart)
end

@testitem "divergence" begin
    using CompScienceMeshes, LinearAlgebra

    T = Float64
    s = simplex(
        point(3,0,0),
        point(0,2,1),
        point(-1,-1,-1),
    )

    order = 2
    ϕ = BEAST.GWPDivRefSpace{T,order}()
    
    u = T(0.2432); dx = sqrt(eps(T))
    v = T(0.5786); dy = sqrt(eps(T))

    # p00 = neighborhood(s, (u,v))
    # p10 = neighborhood(s, (u+du,v))
    # p01 = neighborhood(s, (u, v+dv))
    
    p00 = neighborhood(s, (u,v))
    t1 = normalize(tangents(p00,1))
    t2 = normal(p00) × t1

    p10 = neighborhood(s, carttobary(s, cartesian(p00) + dx*t1 + 0*t2))
    p01 = neighborhood(s, carttobary(s, cartesian(p00) + 0*t1 + dy*t2))

    ϕ00 = ϕ(p00)
    ϕ10 = ϕ(p10)
    ϕ01 = ϕ(p01)


    for (f00, f10, f01) in zip(ϕ00, ϕ10, ϕ01)
        df1dx = dot(t1, f10.value - f00.value) / dx
        df2dy = dot(t2, f01.value - f00.value) / dy
        curl_num = df1dx + df2dy
        curl_ana = f00.divergence
        # @show curl_num,  curl_ana, abs(curl_num - curl_ana)
        @test curl_num ≈ curl_ana atol=sqrt(eps(T))*100
    end
end