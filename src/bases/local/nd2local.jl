"""
Local shape function `r` is the one whose field lines straddle local edge `r`,
which is the edge adjacent to vertex `r`.

This is not the edge starting at vertex `r`. The downside of this local numbering
scheme is that it cannot be extended to cells that are not simplices because
there is no well defined concept of adjacent-ness.
"""
mutable struct ND2RefSpace{T} <: RefSpace{T,8} end

function (ϕ::ND2RefSpace)(nbd)

    u, v = parametric(nbd)
    n = normal(nbd)
    j = jacobian(nbd)

    tu = tangents(nbd,1)
    tv = tangents(nbd,2)

    d = 2/j

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
        (value = n × -(-8*b7-8*b8+12*b3+6*b5-4*b1+6*b6)*inv_j, curl = (-24*u-24*v+18)*inv_j),
        (value = n × -(8*b8-2*b3-6*b5+2*b1-4*b6)*inv_j, curl = (8*u+16*v-6)*inv_j),
        (value = n × (4*b3-8*b7+6*b4+2*b6-2*b2)*inv_j, curl = (6-24*u)*inv_j),
        (value = n × (8*b7+8*b8-6*b3-6*b4-12*b6+4*b2)*inv_j, curl = (24*u+24*v-18)*inv_j),
        (value = n × (2*b3-8*b8+4*b6)*inv_j, curl = (6-24*v)),
        (value = n × (4*b3-8*b7+2*b6)*inv_j, curl = (6-24*u)*inv_j),
        (value = n × (-16*b7-8*b8+16*b3+8*b6)*inv_j, curl = (-48*u-24*v+24)*inv_j),
        (value = n × (-8*b7-16*b8+8*b3+16*b6)*inv_j, curl = (-32*u-48*v+24)*inv_j)
    ))

end


#= function restrict(ϕ::NDRefSpace{T}, dom1, dom2) where T

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

        # find the tangent in this point to the edge
        t = b - a
        # l = norm(t)
        # n = dom2.normals[1]
        # m = cross(t, n) / l

        # dom1 : the smaller domain on which the reference functions are def'd
        # dom2 : the larger domain on which the function to be restrcited is def'd

        u = carttobary(dom1, c)
        x = neighborhood(dom1, u)

        y = ϕ(x)

        for j in 1:K
            # Q[j,i] = dot(y[j][1], m) * l
            Q[j,i] = dot(y[j][1], t)
        end
    end

    return Q
end =#
