function trace(edge_ch, face_ch, fields)
    T = coordtype(edge_ch)
    atol = sqrt(eps(T))

    face_dom = domain(face_ch)
    edge_verts = nodes(edge_ch)
    face_verts = nodes(face_ch)

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






