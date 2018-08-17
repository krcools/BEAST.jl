struct BDMBasis{T,M,P} <: Space{T}
    geo::M
    fns::Vector{Vector{Shape{T}}}
    pos::Vector{P}
end

refspace(s::BDMBasis{T}) where {T} = BDMRefSpace{T}()
subset(s::BDMBasis,I) = BDMBasis(s.geo, s.fns[I], s.pos[I])

function brezzidouglasmarini(mesh)
    edges = skeleton(mesh, 1)
    cps = cellpairs(mesh, edges, dropjunctionpair=true)
    ids = findall(x -> x>0, cps[2,:])
    brezzidouglasmarini(mesh, cps[:,ids])
end

function brezzidouglasmarini(mesh, cellpairs::Array{Int,2})

    @warn "brezzidouglasmarini(mesh, cellpairs) assumes mesh is oriented"

    @assert size(cellpairs,1) == 2

    T = coordtype(mesh)
    P = vertextype(mesh)

    S = Shape{T}
    F = Vector{Shape{T}}

    nf = 2*size(cellpairs,2)
    fns = Vector{F}(undef, nf)
    pos = Vector{P}(undef, nf)

    for i in axes(cellpairs)[2]
        c1, c2 = cellpairs[:,i]
        cell1 = cells(mesh)[c1]
        cell2 = cells(mesh)[c2]
        e1, e2 = getcommonedge(cell1, cell2)
        @assert e1*e2 < 0
        e1, e2 = abs(e1), abs(e2)
        fns[2*(i-1)+1] = [ S(c1, 2*(e1-1)+1 ,+1.0), S(c2, 2*(e2-1)+2,-1.0)]
        fns[2*(i-1)+2] = [ S(c1, 2*(e1-1)+2 ,+1.0), S(c2, 2*(e2-1)+1,-1.0)]

        v1 = cell1[mod1(e1+1,3)]
        v2 = cell1[mod1(e1+2,3)]
        edge = simplex(mesh.vertices[[v1,v2]]...)
        cntr = cartesian(center(edge))
        pos[2*(i-1)+1] = cntr
        pos[2*(i-1)+2] = cntr
    end

    BDMBasis(mesh, fns, pos)
end
