struct RT2RefSpace{T} <: RefSpace{T,8} end

<<<<<<< HEAD
# valuetype(ref::RTRefSpace{T}, charttype) where {T} = SVector{3,Tuple{SVector{universedimension(charttype),T},T}}
function valuetype(ref::RT2RefSpace{T}, charttype::Type) where {T}
    SVector{universedimension(charttype),T}
end

function (ϕ::RT2RefSpace)(mp)

    u, v = parametric(mp)
    j = jacobian(mp)

    tu = tangents(mp,1)
    tv = tangents(mp,2)

    inv_j = 1/j
    
    b1 = tu
    b2 = tv
    b3 = u*tu
    b4 = u*tv
    b5 = v*tu
    b6 = v*tv
    b7 = u*b3+u*b6
    b8 = v*b3+v*b6

    return SVector((
        (value=(8*b7+8*b8-12*b3-6*b5+4*b1-6*b6)*inv_j, divergence=(24*u+24*v-18)*inv_j),
        (value=(-8*b8+2*b3+6*b5-2*b1+4*b6)*inv_j, divergence=(6-24*v)*inv_j),
        (value=(4*b3-8*b7+6*b4+2*b6-2*b2)*inv_j, divergence=(6-24*u)*inv_j),
        (value=(8*b7+8*b8-6*b3-6*b4-12*b6+4*b2)*inv_j, divergence=(24*u+24*v-18)*inv_j),
        (value=(2*b3-8*b8+4*b6)*inv_j, divergence=(6-24*v)*inv_j),
        (value=(4*b3-8*b7+2*b6)*inv_j, divergence=(6-24*u)*inv_j),        
        (value=(-16*b7-8*b8+16*b3+8*b6)*inv_j, divergence=(-48*u-24*v+24)*inv_j),
        (value=(-8*b7-16*b8+8*b3+16*b6)*inv_j, divergence=(-24*u-48*v+24)*inv_j)
    ))
end

#divergence(ref::RT2RefSpace, sh, el) = Shape(sh.cellid, 1, sh.coeff/volume(el))

"""
    ntrace(refspace, element, localindex, face)

Compute the normal trace of all local shape functions on `elements` belonging to
`refspace` on `face`. This function returns a matrix expressing the traces of local
shape functions in `refspace` as linear combinations of functions in the local
trace space. Cf. `restrict`. `localindex` is the index of `face` in the enumeration
of faces of `elements`. In many special cases knowing this index allows for highly
optimised implementations.
"""
function ntrace(x::RT2RefSpace, el, q, fc)
    t = zeros(scalartype(x),1,3)
    t[q] = 1 / volume(fc)
    return t
end

#= function restrict(ϕ::RT2RefSpace{T}, dom1, dom2) where T

    K = numfunctions(ϕ)
    D = dimension(dom1)

    @assert K == 3
    @assert D == 2
    @assert D == dimension(dom2)

    Q = zeros(T,K,K)
    for i in 1:K

        # find the center of edge i of dom2
        a = dom2.vertices[mod1(i+1,D+1)]
        b = dom2.vertices[mod1(i+2,D+1)]
        c = (a + b) / 2

        # find the outer binormal there
        t = b - a
        l = norm(t)
        n = dom2.normals[1]
        m = cross(t, n) / l

        u = carttobary(dom1, c)
        x = neighborhood(dom1, u)

        y = ϕ(x)

        for j in 1:K
            Q[j,i] = dot(y[j][1], m) * l
        end
    end

    return Q
end =#


const _vert_perms_rt2 = [
    (1,2,3),
    (2,3,1),
    (3,1,2),
    (2,1,3),
    (1,3,2),
    (3,2,1),
]
const _dof_perms_rt2 = [
    (1,2,3,4,5,6,7,8),
    (5,6,1,2,3,4,7,8),
    (3,4,5,6,1,2,7,8),
    (4,3,2,1,6,5,7,8),
    (2,1,6,5,4,3,7,8),
    (6,5,4,3,2,1,7,8),
]

function dof_permutation(::RT2RefSpace, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_rt2)
    @assert i != nothing
    return _dof_perms_rt2[i]
end

function dof_perm_matrix(::RT2RefSpace, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_rt2)
=======
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
        (value=(8*u*v-2*u-6*v+2)*tu + (4*v*(2*v-1))*tv, divergence=inv_j),
        (value=(-8*u^2-8*u*v+12*u+6*v-4)*tu + (2*v*(-4*u-4*v+3))*tv, divergence=inv_j),
        (value=-(2*u*(4*u+4*v-3))*tu - (8*u*v-6*u+8*v^2-12*v+4)*tv, divergence=inv_j),
        (value=-(4*u*(1-2*u))*tu - (-8*u*v+6*u+2*v-2)*tv, divergence=inv_j),
        (value=-(4u*(1-2u))*tu - (2v*(1-4u))*tv, divergence=inv_j),
        (value=-(2u*(1-4v))*tu - (4v*(1-2v))*tv, divergence=inv_j),
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
            # x = cartesian(p_edge)
            # u = carttobary(chart, x)

            p_refedge = neighborhood(refedge,s)
            u = cartesian(p_refedge)
            p_refchart = neighborhood(refchart, u)
            t_refedge = tangents(p_refedge,1)
            m_refedge = point(-t_refedge[2], t_refedge[1])
            m_refedge /= norm(m_refedge)
            q0ref = s[1] * m_refedge
            q1ref = (1-s[1]) * m_refedge

            nxq0ref = point(-q0ref[2], q0ref[1])
            nxq1ref = point(-q1ref[2], q1ref[1])
            p_chart = neighborhood(chart, u)
            n_chart = normal(p_chart)
            J_chart = jacobian(p_chart)
            t1 = tangents(p_chart,1)
            t2 = tangents(p_chart,2)
            q0 = -n_chart × (nxq0ref[1]*t1 + nxq0ref[2]*t2) / J_chart^2
            q1 = -n_chart × (nxq1ref[1]*t1 + nxq1ref[2]*t2) / J_chart^2
            
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

        q6 = -n_chart × (nxq6ref[1] * t1 + nxq6ref[2] * t2) / J_chart^2
        q7 = -n_chart × (nxq7ref[1] * t1 + nxq7ref[2] * t2) / J_chart^2

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
>>>>>>> fce8340c8eb82b7d423926c93897058123a66311
    @assert i != nothing
    return _dof_rt2perm_matrix[i]
end