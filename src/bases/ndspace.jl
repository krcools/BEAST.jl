struct NDBasis{T,M,P} <: Space{T}
    geo::M
    fns::Vector{Vector{Shape{T}}}
    pos::Vector{P}
end

NDBasis(geo, fns) = NDBasis(geo, fns, Vector{vertextype(geo)}(undef,length(fns)))

refspace(s::NDBasis) = NDRefSpace{scalartype(s)}()


function nedelec(surface, edges=skeleton(surface,1))

    T = coordtype(surface)
    # P = eltype(surface.vertices)
    P = vertextype(surface)
    num_edges = numcells(edges)

    C = connectivity(edges, surface, identity)
    rows = rowvals(C)
    vals = nonzeros(C)

    fns = Vector{Vector{Shape{T}}}(undef,num_edges)
    pos = Vector{P}(undef,num_edges)
    for (i,edge) in enumerate(edges)

        fns[i] = Vector{Shape{T}}()
        pos[i] = cartesian(center(chart(edges,edge)))

        for k in nzrange(C,i)

            j = rows[k] # j is the index of a cell adjacent to edge
            s = vals[k] # s contains the oriented (signed) local index of edge[i] in cell[j]

            # i == 3 && @show s
            push!(fns[i], Shape{T}(j, abs(s), sign(s)))
        end
    end

    NDBasis(surface, fns, pos)
end

# function LinearAlgebra.cross(::NormalVector, s::NDBasis)
#     # @assert CompScienceMeshes.isoriented(s.geo)
#     RTBasis(s.geo, s.fns, s.pos)
# end

function LinearAlgebra.cross(::NormalVector, s::NDBasis)
    @assert CompScienceMeshes.isoriented(s.geo)
    fns = similar(s.fns)
    for (i,fn) in pairs(s.fns)
        fns[i] = [Shape(sh.cellid, sh.refid, -sh.coeff) for sh in fn]
    end
    RTBasis(s.geo, fns, s.pos)
end
function LinearAlgebra.cross(s::NDBasis, ::NormalVector)
    @assert CompScienceMeshes.isoriented(s.geo)
    fns = similar(s.fns)
    for (i,fn) in pairs(s.fns)
        fns[i] = [Shape(sh.cellid, sh.refid, sh.coeff) for sh in fn]
    end
    RTBasis(s.geo, fns, s.pos)
end
function curl(space::NDBasis)
    divergence(n × (n × (n × space)))
end
