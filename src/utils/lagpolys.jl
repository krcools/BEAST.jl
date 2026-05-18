function _lagpoly(nodes, i, s, i0=1, i1=length(nodes))
    T = eltype(s)
    r = one(T)
    si = nodes[i]
    for j in i0:i1
        j == i && continue
        sj = nodes[j]
        r *= (s - sj) / (si - sj)
    end
    return r
end

function _lagpoly_diff(nodes, i, s, i0=1, i1=length(nodes))
    T = typeof(s)
    r = zero(T)
    si = nodes[i]
    for p in i0:i1
        p == i && continue
        sp = nodes[p]
        rp = one(T)
        for j in i0:i1
            j == i && continue
            j == p && continue
            sj = nodes[j]
            rp *= (s - sj) / (si - sj)
        end
        rp *= 1 / (si - sp)
        r += rp
    end
    return r
end


# Versions of the above functions to be used with @generated functions. These return expressions 
# that can be evaluated at compile time, which allows for more efficient code when the number of
# nodes is known at compile time.
function gen_lagpoly(nodes, i, s, i0, i1, T)
    p = :(one($T))
    for j in i0:i1
        j == i && continue
        p = :(($p * ($s - $nodes[$j]) / ($nodes[$i] - $nodes[$j])))
    end
    return p
end

function gen_lagpoly_diff(nodes, i, s, i0, i1, T)
    dp = :(zero($T))
    for j in i0:i1
        j == i && continue
        p = :(one($T))
        for k in i0:i1
            (k == i || k == j) && continue
            p = :(($p * ($s - $nodes[$k]) / ($nodes[$i] - $nodes[$k])))
        end
        dp = :(($dp + $p / ($nodes[$i] - $nodes[$j])))
    end
    return dp
end

function gen_sylpoly(nodes, i, s, T)
    gen_lagpoly(nodes, i, s, 1, i, T)
end

function gen_sylpoly_diff(nodes, i, s, T)
    gen_lagpoly_diff(nodes, i, s, 1, i, T)
end

function gen_sylpoly_shift(nodes, i, s, T)
    gen_lagpoly(nodes, i, s, 2, i, T)
end

function gen_sylpoly_shift_diff(nodes, i, s, T)
    gen_lagpoly_diff(nodes, i, s, 2, i, T)
end