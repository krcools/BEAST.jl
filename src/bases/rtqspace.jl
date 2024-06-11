struct RTQSpace{T,M,P} <: Space{T}
    geo::M
    fns::Vector{Vector{Shape{T}}}
    pos::Vector{P}
end

# RTQSpace(g::M, fns::Vector{})

function positions(s::RTQSpace) s.pos end 
function refspace(s::RTQSpace{T}) where {T} RTQuadRefSpace{T}() end
function subset(rt::RTQSpace,I) RTQSpace(rt.geo, rt.fns[I], rt.pos[I]) end

function raviartthomas(
    mesh::CompScienceMeshes.QuadMesh{T},
    edges::CompScienceMeshes.AbstractMesh{3,2,T},
    connectivity, orientations) where {T<:Any}

    fns = Vector{Vector{Shape{T}}}(undef, length(edges))
    pos = Vector{SVector{3,T}}(undef, length(edges))

    rows = rowvals(connectivity)
    vals = nonzeros(connectivity)
    for (i,edge) in enumerate(edges)
        σ = orientations[i]
        fn = map(zip(nzrange(connectivity, i),(σ,-σ))) do (j, α)
            Shape{T}(j, abs(vals[j]), α)
        end
        fns[i] = fn
        pos[i] = cartesian(CompScienceMeshes.center(chart(edges, edge)))
    end

    return RTQSpace(mesh, fns, pos)
end

@testitem "RTQSpace construction" begin
    using CompScienceMeshes

    m = CompScienceMeshes.meshrectangle(2.0, 2.0, 1.0; structured=:quadrilateral)
    edges = skeleton(m, 1)
    edges_int = submesh(!in(boundary(edges)),  edges)

    c = CompScienceMeshes.connectivity(edges_int, m)
    @test size(c) == (length(m), length(edges_int))
    o = ones(length(edges_int))
    s = raviartthomas(m, edges_int, c, o)
    @test numfunctions(s) == 4
end