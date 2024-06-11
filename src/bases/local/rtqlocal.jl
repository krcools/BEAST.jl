struct RTQuadRefSpace{T} <: DivRefSpace{T,3} end

function (ϕ::RTQuadRefSpace{T})(p) where {T}

    u, v = parametric(p)
    j = jacobian(p)

    D = tangents(p)

    eu = point(T,1,0)
    ev = point(T,0,1)

    i = one(T)
    vals = SVector(
        (value=(v-1)*ev, divergence=i),
        (value=u*eu,     divergence=i),
        (value=v*ev,     divergence=i),
        (value=(u-1)*eu, divergence=i),
    )

    map(vals) do f
        (value=D*f.value/j, divergence=f.divergence/j) end
end


function interpolate(fields, interpolant::RTQuadRefSpace{T}, chart) where {T}

    refchart = domain(chart)
    Q = map(zip(faces(chart), faces(refchart))) do (edge,refedge)
        s = T(0.5)
        p_edge = neighborhood(edge, s)
        p_refedge = neighborhood(refedge, s)

        p_refchart = neighborhood(refchart, cartesian(p_refedge))
        p_chart = neighborhood(chart, cartesian(p_refchart))
        n_chart = normal(p_chart)
        
        t_edge = tangents(p_edge, 1)
        m_edge = -cross(t_edge, n_chart)
        
        fieldvals = fields(p_chart)
        map(fv -> dot(fv,m_edge), fieldvals)
    end

    return hcat(Q...)
end

@testitem "interpolate" begin
    using CompScienceMeshes

    p1 = point(0,0,0)
    p2 = point(2,0,0)
    p3 = point(2,4,0)
    p4 = point(0,4,0)

    quad = CompScienceMeshes.Quadrilateral(p1,p2,p3,p4)
    f(p) = (x = cartesian(p); return [point(x[1]+2, -x[2]-3, 0)])

    rtq = BEAST.RTQuadRefSpace{Float64}()
    Q = BEAST.interpolate(f, rtq, quad)
    p = neighborhood(quad, point(0.5, 0.5))
    val1 = f(p)[1]
    val2 = sum(Q[1,i] * ϕ.value for (i,ϕ) in zip(axes(Q,2), rtq(p)))

    @show val1
    @show val2
    @test val1 ≈ val2
end