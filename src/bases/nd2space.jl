struct ND2Basis{T,M,P} <: Space{T}
    geo::M
    fns::Vector{Vector{Shape{T}}}
    pos::Vector{P}
end

ND2Basis(geo, fns) = ND2Basis(geo, fns, Vector{vertextype(geo)}(undef,length(fns)))

refspace(s::ND2Basis) = ND2RefSpace{scalartype(s)}()


function nedelec2(surface, edges=skeleton(surface,1))

    T = coordtype(surface)
    # P = eltype(surface.vertices)
    P = vertextype(surface)
    num_edges = numcells(edges)

    C = connectivity(edges, surface, identity)
    rows = rowvals(C)
    vals = nonzeros(C)

    Cells = cells(surface)
    num_cells = size(Cells,1)

    fns = Vector{Vector{Shape{T}}}(undef,2*(num_edges+num_cells))
    pos = Vector{P}(undef,2*(num_edges+num_cells))

    chooseref1(x) = if x<0 return 0 else return -1 end
    chooseref2(x) = if x<0 return -1 else return 0 end
    for (i,edge) in enumerate(edges)

        fns[2*i-1] = Vector{Shape{T}}()
        fns[2*i] = Vector{Shape{T}}()
        pos[2*i-1] = cartesian(center(chart(edges,edge)))
        pos[2*i] = cartesian(center(chart(edges,edge)))

        sgn = 1.0

        for k in nzrange(C,i)

            j = rows[k] # j is the index of a cell adjacent to edge
            s = vals[k] # s contains the oriented (signed) local index of edge[i] in cell[j]

            # i == 3 && @show s
            push!(fns[2*i-1], Shape{T}(j, 2*abs(s)+chooseref1(sign(s)), T(sgn)))
            sgn *= -1*sgn
        end

        sgn = 1.0

        for k in nzrange(C,i)

            j = rows[k] # j is the index of a cell adjacent to edge
            s = vals[k] # s contains the oriented (signed) local index of edge[i] in cell[j]

            # i == 3 && @show s
            push!(fns[2*i], Shape{T}(j, 2*abs(s)+chooseref2(sign(s)), T(sgn)))
            sgn*= -sgn
        end
    end

    for (i,cell) in enumerate(Cells)    
        fns[2*num_edges+2*i-1] = Vector{Shape{T}}()
        fns[2*num_edges+2*i] = Vector{Shape{T}}()
        push!(fns[2*num_edges+2*i-1], Shape{T}(i, 7, T(1.0)))
        push!(fns[2*num_edges+2*i], Shape{T}(i, 8, T(1.0)))
        ctr = cartesian(center(chart(surface,i)))
        pos[2*num_edges+2*i-1] = ctr
        pos[2*num_edges+2*i] = ctr
    end

    ND2Basis(surface, fns, pos)
end

function LinearAlgebra.cross(::NormalVector, s::ND2Basis)
    # @assert CompScienceMeshes.isoriented(s.geo)
    RT2Basis(s.geo, s.fns, s.pos)
end


function curl(space::ND2Basis)
    divergence(n × space)
end
