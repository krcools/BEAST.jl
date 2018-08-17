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
