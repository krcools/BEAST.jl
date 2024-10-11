struct GWPDivRefSpace{T,Degree} <: RefSpace{T} end

function numfunctions(x::GWPDivRefSpace{<:Any,D},
    dom::CompScienceMeshes.ReferenceSimplex{2}) where {D}
    (D+1)*(D+3)
end
function dimtype(x::GWPDivRefSpace{<:Any,D},
    dom::CompScienceMeshes.ReferenceSimplex{2}) where {D}
    Val{(D+1)*(D+3)}
end

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

@testitem "divergence - pointwise" begin
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


function divergence(localspace::GWPDivRefSpace, sh, ch)
    fns = divergence(localspace, sh.cellid, ch)
    α = sh.coeff
    S = typeof(sh)
    return S[S(s.cellid, s.refid, α*s.coeff) for s in fns[sh.refid]]
end


function divergence(localspace::GWPDivRefSpace, cellid::Int, ch)
    divergence(localspace, cellid, ch, dimtype(localspace, domain(ch)))
end


function divergence(localspace::GWPDivRefSpace{T,D}, cellid::Int, ch, ::Type{Val{N}}) where {N,D,T}
    function fields(p)
        map(localspace(p)) do x
            x.divergence
        end
    end
    atol = sqrt(eps(T))
    Dim = 2
    NFout = div((Dim+1)*(Dim+2), 2)
    lag = LagrangeRefSpace{T,D,Dim+1,NFout}()
    coeffs = interpolate(fields, lag, ch)
    S = BEAST.Shape{T}
    A = Vector{Vector{S}}(undef, size(coeffs,1))
    for r in axes(coeffs,1)
        A[r] = collect(S(cellid, c, coeffs[r,c]) for c in axes(coeffs,2) if abs(c) > atol)
    end
    return SVector{N}(A)
end


@testitem "divergence - chartwise" begin
    using CompScienceMeshes
    const CSM = CompScienceMeshes

    T = Float64
    D = 4
    NF = binomial(2+D, 2)
    gwp = BEAST.GWPDivRefSpace{T,D}()
    lgx = BEAST.LagrangeRefSpace{T,D,3,10}()

    ch = CSM.simplex(
        point(1,0,0),
        point(0,1,0),
        point(0,0,0))

    divfns = BEAST.divergence(gwp, 1, ch)

    p = neighborhood(ch, (0.2, 0.6))
    gwp_vals = gwp(p)
    lgx_vals = lgx(p)

    err = similar(Vector{T}, axes(gwp_vals))
    for i in eachindex(gwp_vals)
        val1 = gwp_vals[i].divergence
        val2 = zero(T)
        for sh in divfns[i]
            val2 += sh.coeff * lgx_vals[sh.refid].value
        end
        err[i] = abs(val1 - val2)
    end
    atol = sqrt(eps(T))
    @test all(err .< atol)
end