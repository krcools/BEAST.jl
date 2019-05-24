"""
Local shape function `r` is the one whose field lines straddle local edge `r`,
which is the edge adjacent to vertex `r`.

This is not the edge starting at vertex `r`. The downside of this local numbering
scheme is that it cannot be extended to cells that are not simplices because
there is no well defined concept of adjacent-ness.
"""
mutable struct NDRefSpace{T} <: RefSpace{T,3} end

function (ϕ::NDRefSpace)(nbd)

    u, v = parametric(nbd)
    n = normal(nbd)
    j = jacobian(nbd)

    tu = tangents(nbd,1)
    tv = tangents(nbd,2)

    d = 2/j

    return SVector((
        (n × (tu*(u-1) + tv*v    ) / j, d),
        (n × (tu*u     + tv*(v-1)) / j, d),
        (n × (tu*u     + tv*v    ) / j, d)
    ))

end


function restrict(ϕ::NDRefSpace{T}, dom1, dom2) where T

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
            Q[i,j] = dot(y[j][1], t)
        end
    end

    return Q
end
