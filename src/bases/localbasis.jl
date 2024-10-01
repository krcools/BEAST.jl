abstract type RefSpace{T} end
abstract type DivRefSpace{T} <: RefSpace{T} end

function pushforwardcurl(vals, nbhd)
    tu = tangents(nbhd, 1)
    tv = tangents(nbhd, 2)
    j = jacobian(nbhd)
    n = normal(nbhd)
    map(vals) do v
        f = v.value
        σ = v.curl
        gu = +f[2]
        gv = -f[1]
        (value=cross(n, (gu*tu + gv*tv)/j), curl=σ/j)
    end
end


function pushforwarddiv(vals, nbhd)
    tu = tangents(nbhd, 1)
    tv = tangents(nbhd, 2)
    j = jacobian(nbhd)
    n = normal(nbhd)
    map(vals) do v
        f = v.value
        σ = v.curl
        (value=(f[1]*tu + f[2]*tv)/j, curl=σ/j)
    end
end