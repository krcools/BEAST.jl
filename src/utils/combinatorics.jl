# given a simplex and a face returns:
# +v if face is the v-th face of the simplex oriented according to the simplex
# -v if face is the v-th face of the simplex oriented oppositely to the simplex
# 0 is face is not a face of the simplex
function relorientation(face, simplex)

    v = setdiff(simplex, face)
    length(v) == 1 || return 0

    # find the position of the missing vertex
    v = v[1]
    i = Base.something(findfirst(isequal(v),simplex),0)
    s = (-1)^(i-1)

    # remove that vertex from the simplex
    face2 = Array{Int}(undef,length(simplex)-1)
    for j in 1 : i-1
        face2[j] = simplex[j]
    end
    for j in i : length(simplex)-1
        face2[j] = simplex[j+1]
    end

    # get the permutation that maps face to face2
    #p = indexin(face, face2)
    p = [ something(findfirst(isequal(v),face2),0) for v in face ]

    return s * levicivita(p) * i
end


"""
    getcommonedge(cell1, cell2) -> e1, e2, edge

Returns in edge the common vertices of cell1 and cell2. e1 contains the index
of the vertex of cell1 opposite to this common edge, and with a plus or minus
sign depending on whether the orientation of the common edge is along or
against the internal orientation of cell1. Similar for e2.
"""
function getcommonedge(cell1, cell2)
    isct = intersect(cell1, cell2)
    relorientation(isct, cell1), relorientation(isct, cell2), isct
end
