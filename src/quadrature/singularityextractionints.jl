abstract type SingularityExtractionRule end
regularpart_quadrule(qr::SingularityExtractionRule) = qr.regularpart_quadrule

function momintegrals!(op,
    g, f,t, s,
    z, qrule::SingularityExtractionRule)

    womps = qrule.outer_quad_points

    sop = singularpart(op)
    rop = regularpart(op)

    regqrule = regularpart_quadrule(qrule)
    momintegrals!(rop, g, f, t, s, z, regqrule)

    for p in 1 : length(womps)
        x = womps[p].point
        dx = womps[p].weight

        innerintegrals!(sop, x, g, f, t, s, z, qrule, dx)
    end
end
