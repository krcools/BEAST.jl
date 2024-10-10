struct RTQuadRefSpace{T} <: DivRefSpace{T} end

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

function numfunctions(ϕ::RTQuadRefSpace, dom::CompScienceMeshes.RefQuadrilateral) 4 end

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


const _vrtperm_matrix_rtq = [
    (1,2,3,4),
    (2,3,4,1),
    (3,4,1,2),
    (4,1,2,3),
    (2,1,4,3),
    (1,4,3,2),
    (4,3,2,1),
    (3,2,1,4),
]

const _dofperm_matrix_rtq = [
    @SMatrix[1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0],
    @SMatrix[0.0 0.0 0.0 1.0; 1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0],
    @SMatrix[0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0; 1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0],
    @SMatrix[0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0; 1.0 0.0 0.0 0.0],
    @SMatrix[1.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0; 0.0 0.0 1.0 0.0; 0.0 1.0 0.0 0.0],
    @SMatrix[0.0 0.0 0.0 1.0; 0.0 0.0 1.0 0.0; 0.0 1.0 0.0 0.0; 1.0 0.0 0.0 0.0],
    @SMatrix[0.0 0.0 1.0 0.0; 0.0 1.0 0.0 0.0; 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0],
    @SMatrix[0.0 1.0 0.0 0.0; 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0; 0.0 0.0 1.0 0.0],
]

function dof_perm_matrix(::RTQuadRefSpace, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vrtperm_matrix_rtq)
    @assert i != nothing
    return _dofperm_matrix_rtq[i]
end

@testitem "restrict RTQ0" begin
    using CompScienceMeshes
    using Combinatorics

    ref_vertices = [
        point(0,0),
        point(1,0),
        point(1,1),
        point(0,1)]
    vertices = [
        point(0,0,0),
        point(1,0,0),
        point(1,1,0),
        point(0,1,0)]
    chart1 = CompScienceMeshes.Quadrilateral(vertices...)
    for I in BEAST._vrtperm_matrix_rtq
        @show I
        chart2 = CompScienceMeshes.Quadrilateral(
            vertices[I[1]],
            vertices[I[2]],
            vertices[I[3]],
            vertices[I[4]])
        chart2tochart1 = CompScienceMeshes.Quadrilateral(ref_vertices[collect(I)]...)
        rs = BEAST.RTQuadRefSpace{Float64}()
        Q1 = BEAST.dof_perm_matrix(rs, I)
        Q2 = BEAST.restrict(rs, chart1, chart2, chart2tochart1)
        @test Q1 ≈ Q1
    end
end