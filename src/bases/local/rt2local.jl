struct RT2RefSpace{T} <: RefSpace{T} end

function (f::RT2RefSpace)(p)

    u, v = parametric(p)
    j = jacobian(p)

    tu = tangents(p,1)
    tv = tangents(p,2)

    inv_j = 1/j

    # return SVector(
    #     (value=-(4u*(1-2u))*tu - (2v*(1-4u))*tv, divergence=inv_j),
    #     (value=-(2u*(1-4v))*tu - (4v*(1-2v))*tv, divergence=inv_j),
    #     (value=(-8*u^2-8*u*v+12*u+6*v-4)*tu + (2*v*(-4*u-4*v+3))*tv, divergence=inv_j),
    #     (value=(8*u*v-2*u-6*v+2)*tu + (4*v*(2*v-1))*tv, divergence=inv_j),
    #     (value=-(2*u*(4*u+4*v-3))*tu - (8*u*v-6*u+8*v^2-12*v+4)*tv, divergence=inv_j),
    #     (value=-(4*u*(1-2*u))*tu - (-8*u*v+6*u+2*v-2)*tv, divergence=inv_j),
    #     (value=(8*u*(-2*u-v+2))*tu + (8*v*(-2*u-v+1))*tv, divergence=inv_j),
    #     (value=(8*u*(-u-2*v+1))*tu + (8*v*(-u-2*v+2))*tv, divergence=inv_j),
    # )

    return SVector(
        (value=((8*u^2+8*u*v-12*u-6*v+4)*tu + (2*v*(4*u+4*v-3))*tv)*inv_j, divergence=(24*u+24*v-18)*inv_j),
        (value=((-8*u*v+2*u+6*v-2)*tu + (4*v*(-2*v+1))*tv)*inv_j, divergence=(6-24*v)*inv_j),
        (value=((4*u*(1-2*u))*tu + (-8*u*v+6*u+2*v-2)*tv)*inv_j, divergence=(6-24*u)*inv_j),
        (value=((2*u*(4*u+4*v-3))*tu + (8*u*v-6*u+8*v^2-12*v+4)*tv)*inv_j, divergence=(24*u+24*v-18)*inv_j),
        (value=((2u*(1-4v))*tu + (4v*(1-2v))*tv)*inv_j, divergence=(6-24*v)*inv_j),
        (value=((4u*(1-2u))*tu + (2v*(1-4u))*tv)*inv_j, divergence=(6-24*u)*inv_j),
        (value=((8*u*(-2*u-v+2))*tu + (8*v*(-2*u-v+1))*tv)*inv_j, divergence=(-48*u-24*v+24)*inv_j),
        (value=((8*u*(-u-2*v+1))*tu + (8*v*(-u-2*v+2))*tv)*inv_j, divergence=(-24*u-48*v+24)*inv_j)
    )
end

numfunctions(x::RT2RefSpace, dom::CompScienceMeshes.ReferenceSimplex{2}) = 8

function interpolate(fields, interpolant::BEAST.RT2RefSpace, chart)

    T = coordtype(chart)

    Q = Any[]
    refchart = CompScienceMeshes.domain(chart)
    nfields = length(fields(center(chart)))

    for (edge, refedge) in zip(faces(chart), faces(refchart))
        l0 = zeros(T,nfields)
        l1 = zeros(T,nfields)
        qps = CompScienceMeshes.quadpoints(edge,4)
        for (p_edge,w) in qps
            s = parametric(p_edge)

            p_refedge = neighborhood(refedge,s)
            u = cartesian(p_refedge)
            p_refchart = neighborhood(refchart, u)
            t_refedge = tangents(p_refedge,1)
            m_refedge = point(-t_refedge[2], t_refedge[1])
            m_refedge /= norm(m_refedge)
            q0ref = (-1) * (1-s[1]) * m_refedge
            q1ref = (-1) * s[1] * m_refedge

            nxq0ref = point(-q0ref[2], q0ref[1])
            nxq1ref = point(-q1ref[2], q1ref[1])
            p_chart = neighborhood(chart, u)
            n_chart = normal(p_chart)
            J_chart = jacobian(p_chart)
            t1 = tangents(p_chart,1)
            t2 = tangents(p_chart,2)
            q0 = -n_chart × (nxq0ref[1]*t1 + nxq0ref[2]*t2) / J_chart
            q1 = -n_chart × (nxq1ref[1]*t1 + nxq1ref[2]*t2) / J_chart
            
            vals = fields(p_chart)
            J_edge = jacobian(p_edge)
            J_refedge = jacobian(p_refedge)
            l0 .+= [w * dot(f,q0) * J_chart / J_edge * J_refedge for f in vals]
            l1 .+= [w * dot(f,q1) * J_chart / J_edge * J_refedge for f in vals]
        end
        push!(Q,l0)
        push!(Q,l1)
    end

    l6 = zeros(T,nfields)
    l7 = zeros(T,nfields)
    qps = CompScienceMeshes.quadpoints(chart, 4)
    for (p,w) in qps
        
        q6ref = point(1,0)
        q7ref = point(0,1)

        nxq6ref = point(-q6ref[2], q6ref[1])
        nxq7ref = point(-q7ref[2], q7ref[1])

        J_chart = jacobian(p)
        n_chart = normal(p)

        t1 = tangents(p,1)
        t2 = tangents(p,2)

        q6 = -n_chart × (nxq6ref[1] * t1 + nxq6ref[2] * t2) / J_chart
        q7 = -n_chart × (nxq7ref[1] * t1 + nxq7ref[2] * t2) / J_chart

        vals = fields(p)
        l6 .+= [w * dot(f,q6) for f in vals]
        l7 .+= [w * dot(f,q7) for f in vals]
    end

    push!(Q,l6)
    push!(Q,l7)

    return hcat(Q...)
end


function dof_perm_matrix(::RT2RefSpace, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_rt)
    @assert i != nothing
    return _dof_rt2perm_matrix[i]
end

const _dof_rt2perm_matrix = [
    @SMatrix[1 0 0 0 0 0 0 0;  #1. {1,2,3}
                0 1 0 0 0 0 0 0;
                0 0 1 0 0 0 0 0;
                0 0 0 1 0 0 0 0;
                0 0 0 0 1 0 0 0;
                0 0 0 0 0 1 0 0;
                0 0 0 0 0 0 1 0;
                0 0 0 0 0 0 0 1],

    @SMatrix[0 0 0 0 1 0 0 0;  #2. {2,3,1}
                0 0 0 0 0 1 0 0;
                1 0 0 0 0 0 0 0;
                0 1 0 0 0 0 0 0;
                0 0 1 0 0 0 0 0;
                0 0 0 1 0 0 0 0;
                0 0 0 0 0 0 0 -1;
                0 0 0 0 0 0 1 -1],

    @SMatrix[0 0 1 0 0 0 0 0;  #3. {3,1,2}
                0 0 0 1 0 0 0 0;
                0 0 0 0 1 0 0 0;
                0 0 0 0 0 1 0 0;
                1 0 0 0 0 0 0 0;
                0 1 0 0 0 0 0 0;
                0 0 0 0 0 0 -1 1;
                0 0 0 0 0 0 -1 0],

    @SMatrix[0 0 0 1 0 0 0 0;  #4. {2,1,3}
                0 0 1 0 0 0 0 0;
                0 1 0 0 0 0 0 0;
                1 0 0 0 0 0 0 0;
                0 0 0 0 0 1 0 0;
                0 0 0 0 1 0 0 0;
                0 0 0 0 0 0 0 1;
                0 0 0 0 0 0 1 0],

    @SMatrix[0 1 0 0 0 0 0 0;  #5. {1,3,2}
                1 0 0 0 0 0 0 0;
                0 0 0 0 0 1 0 0;
                0 0 0 0 1 0 0 0;
                0 0 0 1 0 0 0 0;
                0 0 1 0 0 0 0 0;
                0 0 0 0 0 0 1 -1;
                0 0 0 0 0 0 0 -1],

    @SMatrix[0 0 0 0 0 1 0 0;  #6. {3,2,1}
                0 0 0 0 1 0 0 0;
                0 0 0 1 0 0 0 0;
                0 0 1 0 0 0 0 0;
                0 1 0 0 0 0 0 0;
                1 0 0 0 0 0 0 0;
                0 0 0 0 0 0 -1 0;
                0 0 0 0 0 0 -1 1]
]

@testitem "restrict RT2" begin
    using CompScienceMeshes
    using LinearAlgebra
    
    ref_vertices = [
        point(1,0),
        point(0,1),
        point(0,0),
    ]
    vertices = [
        point(1,0,0),
        point(0,1,0),
        point(0,0,0),
    ]
    chart1 = simplex(vertices...)
    for I in BEAST._vert_perms_rt
        chart2 = simplex(
            chart1.vertices[I[1]],
            chart1.vertices[I[2]],
            chart1.vertices[I[3]],)
        chart2tochart1 = CompScienceMeshes.simplex(ref_vertices[collect(I)]...)
        rs = BEAST.RT2RefSpace{Float64}()
        Q1 = BEAST.dof_perm_matrix(rs, I)
        Q3 = BEAST.interpolate(rs, chart2, rs, chart1, chart2tochart1)
        @test Q1 ≈ Q3
        @show norm(Q1-Q3)
    end
end