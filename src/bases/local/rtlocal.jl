struct RTRefSpace{T} <: RefSpace{T,3} end

valuetype(ref::RTRefSpace{T}, charttype) where {T} = SVector{3,Tuple{Pt{universedimension(charttype),T},T}}

function (ϕ::RTRefSpace)(mp)

    u, v = parametric(mp)
    j = jacobian(mp)

    tu = tangents(mp,1)
    tv = tangents(mp,2)

    d = 2/j

    return SVector((
        (value=(tu*(u-1) + tv*v    )/j, divergence=d),
        (value=(tu*u     + tv*(v-1))/j, divergence=d),
        (value=(tu*u     + tv*v    )/j, divergence=d)
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
end
