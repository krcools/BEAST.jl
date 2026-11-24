function trace_vertex_permutation(edge_verts, face_verts, atol)
    allunique(face_verts) ||
        throw(ArgumentError("the containing cell has coincident vertices"))
    perm = map(edge_verts) do v
        # Exact matches remain distinct even when several vertices of a thin cell
        # lie within atol. No shared mesh or vertex numbering is required.
        exact = findfirst(==(v), face_verts)
        isnothing(exact) || return exact

        # hypot avoids squaring small distances in Float16/Float32.
        distances = [hypot((v - w)...) for w in face_verts]
        distance, i = findmin(distances)
        isfinite(distance) && distance <= atol ||
            throw(ArgumentError("a trace vertex does not match a vertex of its cell"))
        count(==(distance), distances) == 1 ||
            throw(ArgumentError("a trace vertex has equally close matches"))
        return i
    end
    allunique(perm) ||
        throw(ArgumentError("distinct trace vertices map to the same cell vertex"))
    return perm
end

function trace(edge_ch, face_ch, fields)
    T = coordtype(edge_ch)

    face_dom = domain(face_ch)
    edge_verts = nodes(edge_ch)
    face_verts = nodes(face_ch)

    # The vertices of edge_ch are vertices of face_ch, identified by their positions.
    # The tolerance is relative to the size of face_ch, so that the identification
    # does not depend on the unit of length: with the absolute tolerance sqrt(eps),
    # all vertices of a cell smaller than about 1e-8 matched its first vertex.
    diameter = maximum(norm(a - b) for a in face_verts, b in face_verts)
    atol = sqrt(eps(T)) * diameter
    perm = trace_vertex_permutation(edge_verts, face_verts, atol)
    injection = simplex(vertices(face_dom)[perm])

    function f(s)
        u = neighborhood(injection, s)
        p = neighborhood(face_ch, cartesian(u))
        fields(p)
    end

    return f
end






