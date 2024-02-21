struct RT2RefSpace{T} <: RefSpace{T,8} end

function (f::RT2RefSpace)(p)

    u, v = parametric(mp)
    j = jacobian(mp)

    tu = tangents(mp,1)
    tv = tangents(mp,2)

    inv_j = 1/j

    return SVector(
        (value=(4u*(1-2u))*tu + (2v*(1-4u))*tv, divergence=inv_j),
        (value=(2u*(1-4v))*tu + (4v*(1-2v))*tv, divergence=inv_j),
        (value=(-8*u^2-8*u*v+12*u+6*v-4)*tu + (2*v*(-4*u-4*v+3))*tv, divergence=inv_j),
        (value=(8*u*v-2*u-6*v+2)*tu + (4*v*(2*v-1))*tv, divergence=inv_j),
        (value=(2*u*(4*u+4*v-3))*tu + (8*u*v-6*u+8*v^2-12*v+4)*tv, divergence=inv_j),
        (value=(4*u*(1-2*u))*tu + (-8*u*v+6*u+2*v-2)*tv, divergence=inv_j),
        (value=(8*u*(-2*u-v+2))*tu + (8*v*(-2*u-v+1))*tv, divergence=inv_j),
        (value=(8*u*(-u-2*v+1))*tu + (8*v*(-u-2*v+2))*tv, divergence=inv_j),
    )
end


function interpolate(fields, interpolant::BEAST.RT2RefSpace, chart)

    T = coordtype(chart)

    Q = Any[]
    refchart = CompScienceMeshes.domain(chart).simplex
    nfields = length(fields(center(chart)))

    for (edge, refedge) in zip(faces(chart), faces(refchart))
        l0 = zeros(T,nfields)
        l1 = zeros(T,nfields)
        qps = CompScienceMeshes.quadpoints(edge,4)
        for (p_edge,w) in qps
            s = parametric(p_edge)
            x = cartesian(p_edge)
            u = carttobary(x, chart)

            p_refchart = neighborhood(refchart, u)
            t_refedge = tangents(refedge,1)
            m_refedge = point(-t_refedge[2], t_refedge[1])
            m_refedge /= norm(m_refedge)
            q0ref = s * m_refedge
            q1ref = (1-s) * m_refedge

            nxq0ref = point(-q0ref[2], q0ref[1])
            nxq1ref = point(-q1ref[2], q1ref[1])
            p_chart = neighborhood(u, chart)
            n_chart = normal(p_chart)
            J_chart = jacobian(p_chart)
            t1 = tangents(p_chart,1)
            t2 = tangents(p_chart,2)
            q0 = -n_chart × (nxq0ref[1]*t1 + xnq0ref[2]*t2) / J_chart^2
            q1 = -n_chart × (nxq1ref[1]*t1 + xnq1ref[2]*t2) / J_chart^2
            
            F = fields(p_chart)
            p_refedge = neighborhood(s, refedge)
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

        q6 = -n_chart × (nxq6ref[1] * t1 + nxq6ref[2] * t2) / J_chart^2
        q7 = -n_chart × (nxq7ref[1] * t1 + nxq7ref[2] * t2) / J_chart^2

        vals = fields(p)
        l6 .+= [w * dot(f,q6) for f in vals]
        l7 .+= [w * dot(f,q6) for f in vals]
    end

    return hcat(Q...)
end