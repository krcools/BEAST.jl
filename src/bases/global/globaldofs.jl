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
    perm = [something(findfirst(w -> isapprox(v, w; atol), face_verts), 0)
        for v in edge_verts]
    injection = simplex(vertices(face_dom)[perm])

    function f(s)
        u = neighborhood(injection, s)
        p = neighborhood(face_ch, cartesian(u))
        fields(p)
    end

    return f
end






