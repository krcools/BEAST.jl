function Base.:*(A::AbstractArray, s::Space) where {Space}

    @assert ndims(A) == 2
    @assert size(A,2) == numfunctions(s)

    F = eltype(s.fns)
    fns = [similar(F,0) for i in axes(A,1)]
    geo = s.geo

    for i in axes(A,1)
        for j in axes(A,2)
            A[i,j] == 0 && continue
            shapes = [Shape(sh.cellid, sh.refid, A[i,j]*sh.coeff) for sh in s.fns[j]]
            append!(fns[i], shapes)
        end
        fns[i] = collapse_shapes(fns[i])
    end

    pos = fill(zero(eltype(s.pos)), length(fns))
    return Space(geo, fns, pos)
end


function collapse_shapes(shapevec)
    ids = unique([(sh.cellid,sh.refid) for sh in shapevec])

    S = eltype(shapevec)
    collapsed = [S(id..., 0) for id in ids]

    for sh in shapevec
        pos = findfirst(isequal((sh.cellid, sh.refid)), ids)
        @assert pos != nothing
        old = collapsed[pos]
        @assert old.cellid == sh.cellid
        @assert old.refid == sh.refid
        collapsed[pos] = S(old.cellid, old.refid, old.coeff+sh.coeff)
    end

    collapsed
end
