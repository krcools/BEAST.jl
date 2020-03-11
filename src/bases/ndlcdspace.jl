mutable struct NDLCDBasis{T,M,P} <: Space{T}
    geo::M
    fns::Vector{Vector{Shape{T}}}
    pos::Vector{P}
end

NDLCDBasis(geo, fns) = NDLCDBasis(geo, fns, Vector{vertextype(geo)}(undef,length(fns)))

refspace(space::NDLCDBasis{T}) where {T} = NDLCDRefSpace{T}()
Base.similar(space::NDLCDBasis{T,M,P} where {T,M,P}, geo, fns, pos) = NDLCDBasis(geo, fns, pos)

function nedelecd3d(mesh, faces)
    T = coordtype(mesh)
    P = vertextype(mesh)
    num_faces = numcells(faces)

    C = connectivity(faces, mesh, identity)
    rows = rowvals(C)
    vals = nonzeros(C)

    fns = Vector{Vector{Shape{T}}}(undef,num_faces)
    pos = Vector{P}(undef,num_faces)
    for (i,face) in enumerate(cells(faces))

        fns[i] = Vector{Shape{T}}()
        pos[i] = cartesian(center(chart(faces,face)))

        for k in nzrange(C,i)

            j = rows[k]
            s = vals[k]
            push!(fns[i], Shape{T}(j, abs(s), sign(s)))
        end
    end

    NDLCDBasis(mesh, fns, pos)
end

function nedelecd3d(mesh)
    faces = skeleton(mesh,2)
    nedelecd3d(mesh, faces)
end

ntrace(X::NDLCDBasis, geo, fns) = LagrangeBasis{0,-1,1}(geo, fns, deepcopy(X.pos))
ttrace(X::NDLCDBasis, geo, fns) = NDBasis{}(geo, fns, deepcopy(X.pos))

divergence(space::NDLCDBasis, geo, fns) = LagrangeBasis{0,-1,1}(geo, fns, space.pos)
