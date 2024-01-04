struct RT2RefSpace{T} <: RefSpace{T,8} end

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
    d = 2 * inv_j

    b1 = tu
    b2 = tv
    b3 = u*tu
    b4 = u*tv
    b5 = v*tu
    b6 = v*tv
    b7 = u*b3+u*b6
    b8 = v*b3+v*b6

    return SVector((
        (value=-(-8*b7-8*b8+12*b3+6*b5-4*b1+6*b6)*inv_j, divergence=(-24*u-24*v+18)*inv_j),
        (value=-(8*b8-2*b3-6*b5+2*b1-4*b6)*inv_j, divergence=(8*u+16*v-6)*inv_j),
        (value=(4*b3-8*b7+6*b4+2*b6-2*b2)*inv_j, divergence=(6-24*u)*inv_j),
        (value=(8*b7+8*b8-6*b3-6*b4-12*b6+4*b2)*inv_j, divergence=(24*u+24*v-18)*inv_j),
        (value=(2*b3-8*b8+4*b6)*inv_j, divergence=(6-24*v)),
        (value=(4*b3-8*b7+2*b6)*inv_j, divergence=(6-24*u)*inv_j),
        (value=(-16*b7-8*b8+16*b3+8*b6)*inv_j, divergence=(-48*u-24*v+24)*inv_j),
        (value=(-8*b7-16*b8+8*b3+16*b6)*inv_j, divergence=(-32*u-48*v+24)*inv_j)
    ))
end

divergence(ref::RT2RefSpace, sh, el) = Shape(sh.cellid, 1, sh.coeff/volume(el))

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
    (1,2,3,4,5,6,7,8),
    (1,2,3,4,5,6,7,8),
    (1,2,3,4,5,6,7,8),
    (1,2,3,4,5,6,7,8),
    (1,2,3,4,5,6,7,8),
]

function dof_permutation(::RT2RefSpace, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_rt2)
    @assert i != nothing
    return _dof_perms_rt2[i]
end