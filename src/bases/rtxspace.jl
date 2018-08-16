abstract type Continuity{T} end

function raviartthomas(mesh, ::Type{Continuity{:none}})

    @assert dimension(mesh) == 2

    P = vertextype(mesh)
    S = Shape{coordtype(mesh)}
    F = Vector{S}

    nf = 3*numcells(mesh)
    functions = Vector{F}(undef, nf)
    positions = Vector{P}(undef, nf)

    for i in 1:numcells(mesh)
        functions[3*(i-1)+1] = [S(i,1,+1.0)]
        functions[3*(i-1)+2] = [S(i,2,+1.0)]
        functions[3*(i-1)+3] = [S(i,3,+1.0)]

        ctr = sum(mesh.vertices[mesh.faces[i]])/3
        positions[3*(i-1)+1] = ctr
        positions[3*(i-1)+2] = ctr
        positions[3*(i-1)+3] = ctr
    end

    return RTBasis(mesh, functions, positions)
end
