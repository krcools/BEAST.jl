
function DofInterpolate(basis::LagrangeBasis, field) 

    T = promote_type(scalartype(basis), scalartype(field))

    num_bfs = numfunctions(basis)

    res = Vector{T}(undef, num_bfs)

    for  b in 1 : num_bfs
        bfs = basis.fns[b]

        shape = bfs[1]

        cellid = shape.cellid
        refid = shape.refid

        cell = cells(basis.geo)[cellid]

        tria = chart(basis.geo, cellid)

        vid =  cell[refid]

        if refid == 1 
            v = neighborhood(tria, [1, 0])
        elseif refid == 2
            v = neighborhood(tria, [0, 1])
        else 
            v = neighborhood(tria, [0, 0])
        end

        res[b] = field(v)
    end

    return res
end


function DofInterpolate(basis::LagrangeBasis{D,C,M}, field) where {D,C,U,M<:CompScienceMeshes.Mesh{U,2}}

    T = promote_type(scalartype(basis), scalartype(field))

    num_bfs = numfunctions(basis)

    res = Vector{T}(undef, num_bfs)

    for  b in 1 : num_bfs
        bfs = basis.fns[b]

        res[b] = field(basis.pos[b])
    end

    return res
end


### Piecewise constant elements require separate treatment
### TODO: Probably a rewrite is advisable to also take into account
### dual elements properly.
function DofInterpolate(basis::LagrangeBasis{0,-1,M,T,NF,P}, field) where {M, T, NF, P}

    TT = promote_type(scalartype(basis), scalartype(field))

    num_bfs = numfunctions(basis)

    res = Vector{TT}(undef, num_bfs)

    for  b in 1 : num_bfs
        bfs = basis.fns[b]

        basis.pos[b]

        shape = bfs[1]

        cellid = shape.cellid

        tria = chart(basis.geo, cellid)

        v = neighborhood(tria, [1/3, 1/3])

        res[b] = field(v)
    end

    return res
end

function DofInterpolate(basis::RTBasis, field)

    T = promote_type(scalartype(basis), scalartype(field))

    num_bfs = numfunctions(basis)

    res = Vector{T}(undef, num_bfs)

    for b in 1 : num_bfs

        bfs = basis.fns[b]

        shape = bfs[1]

        cellid = shape.cellid
        refid = shape.refid
        coeff = shape.coeff

        cell = cells(basis.geo)[cellid]

        e = refid
        
        v1 = cell[mod1(e+1,3)]
        v2 = cell[mod1(e+2,3)]

        edge = simplex(basis.geo.vertices[[v1,v2]]...)
        t = tangents(center(edge),1)
        tria = chart(basis.geo, cellid)
        
        n = normal(center(tria))

        bn = normalize(cross(n,t))
      
        v = center(edge)

        res[b] = volume(edge)*coeff*dot(field(v),bn)

    end

    return [r for r in res]

end

function DofInterpolate(basis::BDMBasis, field)

    T = promote_type(scalartype(basis), scalartype(field))

    num_bfs = numfunctions(basis)

    res = Vector{T}(undef, num_bfs)

    for b in 1 : num_bfs
        
        bfs = basis.fns[b]

        shape = bfs[2]

        cellid = shape.cellid
        refid = shape.refid
        coeff = shape.coeff

        cell = cells(basis.geo)[cellid]

        e = (refid-1)÷2+1
        
        v1 = cell[mod1(e+1,3)]
        v2 = cell[mod1(e+2,3)]

        edge = simplex(basis.geo.vertices[[v1,v2]]...)
        t = tangents(center(edge),1)
        tria = chart(basis.geo, cellid)
        
        n = normal(center(tria))

        bn = normalize(cross(n,t))

        if refid == 1 
            v = neighborhood(tria, [0, 1])
        elseif refid == 2
            v = neighborhood(tria, [0, 0])
        elseif refid == 3
            v = neighborhood(tria, [0, 0])
        elseif refid == 4
            v = neighborhood(tria, [1, 0])
        elseif refid == 5
            v = neighborhood(tria, [1, 0])
        else
            v = neighborhood(tria, [0, 1])
        end

        res[b] = volume(edge)*coeff*dot(field(v),bn)

    end

    return res

end

function DofInterpolate(basis::NDLCDBasis, field)

    T = promote_type(scalartype(basis), scalartype(field))

    num_bfs = numfunctions(basis)

    res = Vector{T}(undef, num_bfs)

    for b in 1 : num_bfs
        
        bfs = basis.fns[b]

        shape = bfs[1]

        cellid = shape.cellid
        refid = shape.refid
        coeff = shape.coeff

        cell = cells(basis.geo)[cellid]

        f = refid
        
        v1 = cell[mod1(f+1,4)]
        v2 = cell[mod1(f+2,4)]
        v3 = cell[mod1(f+3,4)]

        v4 = cell[mod1(f,4)]
    
        face = simplex(basis.geo.vertices[[v1,v2,v3]]...)

        n = normal(center(face))

        inside = dot(n, basis.geo.vertices[v4]-basis.geo.vertices[v1])

        n *= -sign(inside)
        
        v = center(face)

        res[b] = volume(face)*coeff*dot(field(v),n)

    end

    return res

end

function DofInterpolate(basis::BEAST.BDM3DBasis, field)

    T = promote_type(scalartype(basis), scalartype(field))

    num_bfs = numfunctions(basis)

    res = Vector{T}(undef, num_bfs)

    for b in 1 : num_bfs
        
        bfs = basis.fns[b]

        shape = bfs[1]

        cellid = shape.cellid
        refid = shape.refid
        coeff = shape.coeff

        cell = cells(basis.geo)[cellid]
    
        f = (refid-1)÷3+1
        
        v1 = cell[mod1(f+1,4)]
        v2 = cell[mod1(f+2,4)]
        v3 = cell[mod1(f+3,4)]

        v4 = cell[mod1(f,4)]

        face = simplex(basis.geo.vertices[[v1,v2,v3]]...)

        n = normal(center(face))

        inside = dot(n, basis.geo.vertices[v4]-basis.geo.vertices[v1])

        n *= -sign(inside)

        tet = simplex(basis.geo.vertices[cell]...)

        if refid == 1 
            v = neighborhood(tet, [0, 1, 0])
        elseif refid == 2
            v = neighborhood(tet, [0, 0, 1])
        elseif refid == 3
            v = neighborhood(tet, [0, 0, 0])
        elseif refid == 4
            v = neighborhood(tet, [0, 0, 1])
        elseif refid == 5
            v = neighborhood(tet, [0, 0, 0])
        elseif refid == 6
            v = neighborhood(tet, [1, 0, 0])
        elseif refid == 7
            v = neighborhood(tet, [0, 0, 0])
        elseif refid == 8
            v = neighborhood(tet, [1, 0, 0])
        elseif refid == 9
            v = neighborhood(tet, [0, 1, 0])
        elseif refid == 10
            v = neighborhood(tet, [1, 0, 0])
        elseif refid == 11
            v = neighborhood(tet, [0, 1, 0])
        else
            v = neighborhood(tet, [0, 0, 1])
        end

        res[b] = volume(face)*coeff*dot(field(v),n)

    end

    return res

end

function Centervalues(mesh::Mesh,f::Functional)
    num_cells = numcells(mesh)
    T = coordtype(mesh)
    res = Vector{SVector{3,Complex{T}}}(undef, num_cells)

    for i in 1:num_cells
        cell = cells(mesh)[i]

        tet = chart(mesh, cell)

        v = center(tet)

        res[i] = f(v)

    end

    return res

end

function EvalCenter(basis::Space, coeff)

    mesh = basis.geo
    T = coordtype(mesh)
    num_cells = numcells(mesh)
    num_bfs = numfunctions(basis)
    ref_space = refspace(basis)

    res = Vector{SVector{3,Complex{T}}}(undef, num_cells)

    for i in 1:num_cells
        res[i] = [0,0,0]
    end

    for b in 1 : num_bfs
        
        bfs = basis.fns[b]

        for shape in bfs
       

            cellid = shape.cellid
            refid = shape.refid
            a = shape.coeff

            cell = cells(mesh)[cellid]

            tet = chart(mesh, cell)

            v = center(tet)

            local_bfs = ref_space(v)

            res[cellid] += a*coeff[b]*local_bfs[refid].value
        end

    end

    return res

end
