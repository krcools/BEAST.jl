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
