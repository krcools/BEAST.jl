struct GWPDivSpace{T,M,P} <: Space{T}
    geo::M
    fns::Vector{Vector{Shape{T}}}
    pos::P
    degree::Int
end

function refspace(s::GWPDivSpace{T}) where {T} GWPDivRefSpace{T,s.degree}() end
function subset(s::S,I) where {S<:GWPDivSpace} S(s.geo, s.fns[I], s.pos[I], s.degree) end

function gwpdiv(mesh, edges=nothing; order)

    T = coordtype(mesh)
    S = Shape{T}

    space = BEAST.gwpcurl(mesh, edges; order)
    fns = Vector{Vector{S}}(undef, length(space.fns))
    for (i,fn) in enumerate(space.fns)
        fns[i] = [S(s.cellid, s.refid, -s.coeff) for s in fn]
    end

    return GWPDivSpace(space.geo, fns, space.pos, order)
end

@testitem "GWPcurl global: numfunctions" begin
    using CompScienceMeshes

    h = 0.5
    mesh = meshrectangle(1.0, 1.0, 0.5)
    edges = setminus(skeleton(mesh,1), boundary(mesh))

    order = 2
    gwp = BEAST.gwpdiv(mesh, edges; order=order)

    ne = order+1
    nf = order * (order+1)
    Nt = length(edges)*ne + length(mesh)*nf
    @test numfunctions(gwp) == Nt
end