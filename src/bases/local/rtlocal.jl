struct RTRefSpace{T} <: RefSpace{T,3} end

# valuetype(ref::RTRefSpace{T}, charttype) where {T} = SVector{3,Tuple{SVector{universedimension(charttype),T},T}}
function valuetype(ref::RTRefSpace{T}, charttype::Type) where {T}
    SVector{universedimension(charttype),T}
end

function (ϕ::RTRefSpace)(mp)

    u, v = parametric(mp)
    j = jacobian(mp)

    tu = tangents(mp,1)
    tv = tangents(mp,2)

    inv_j = 1/j
    d = 2 * inv_j

    u_tu = u*tu
    v_tv = v*tv

    return SVector((
        (value=(u_tu-tu + v_tv    )*inv_j, divergence=d),
        (value=(u_tu     + v_tv-tv)*inv_j, divergence=d),
        (value=(u_tu     + v_tv    )*inv_j, divergence=d)
    ))
end

divergence(ref::RTRefSpace, sh, el) = Shape(sh.cellid, 1, sh.coeff/volume(el))

"""
    ntrace(refspace, element, localindex, face)

Compute the normal trace of all local shape functions on `elements` belonging to
`refspace` on `face`. This function returns a matrix expressing the traces of local
shape functions in `refspace` as linear combinations of functions in the local
trace space. Cf. `restrict`. `localindex` is the index of `face` in the enumeration
of faces of `elements`. In many special cases knowing this index allows for highly
optimised implementations.
"""
function ntrace(x::RTRefSpace, el, q, fc)
    t = zeros(scalartype(x),1,3)
    t[q] = 1 / volume(fc)
    return t
end

function restrict(ϕ::RTRefSpace{T}, dom1, dom2) where T

    K = numfunctions(ϕ)
    D = dimension(dom1)

    @assert K == 3
    @assert D == 2
    @assert D == dimension(dom2)

    Q = zeros(T,K,K)
    for i in 1:K

        # find the center of edge i of dom2
        a = verticeslist(dom2)[mod1(i+1,D+1)]
        b = verticeslist(dom2)[mod1(i+2,D+1)]
        c = (a + b) / 2

        # find the outer binormal there
        t = b - a
        l = norm(t)
 #       n = dom2.normals[1]
        n = normalize(tangents(dom2,1)×tangents(dom2,2))
        m = cross(t, n) / l

        u = carttobary(dom1, c)
        x = neighborhood(dom1, u)

        y = ϕ(x)

        for j in 1:K
            Q[j,i] = dot(y[j][1], m) * l
        end
    end

    return Q
end


const _vert_perms_rt = [
    (1,2,3),
    (2,3,1),
    (3,1,2),
    (2,1,3),
    (1,3,2),
    (3,2,1),
]
const _dof_perms_rt = [
    (1,2,3),
    (3,1,2),
    (2,3,1),
    (2,1,3),
    (1,3,2),
    (3,2,1),
]

function dof_permutation(::RTRefSpace, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_rt)
    @assert i != nothing
    return _dof_perms_rt[i]
end