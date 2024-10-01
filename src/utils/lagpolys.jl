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
    r = zero(T)
    si = nodes[i]
    for p in i0:i1
        p == i && continue
        rp = one(T)
        for j in i0:i1
            j == i && continue
            j == p && continue
            sj = nodees[j]
            rp *= (s - sj) / (si - sj)
        end
        rp *= 1 / (si - sp)
        r += rp
    end
    return r
end

function _sylpoly(nodes, i, s)
    _lagpoly(nodes, i, s, 1, i)
end

function _sylpoly_diff(nodes, i, s)
    _lagpoly_diff(nodes, i, s, 1, i)
end

function _sylpoly_shift(nodes, i, s)
    _lagpoly(nodes, i, s, 2, i)
end

function _sylpoly_shift_diff(nodes, i, s)
    _lagpoly_shift(nodes, i, s, 2, i)
end