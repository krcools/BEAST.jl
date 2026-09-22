abstract type SingularityExtractionRule end
regularpart_quadrule(qr::SingularityExtractionRule) = qr.regularpart_quadrule

function integrate!(op,
    g, f,t, s,
    z, qrule::SingularityExtractionRule)

    return integrate!(op, g, f, t, s, z, qrule, quadraturebuffer(qrule, g, f))
end

function integrate!(op,
    g, f,t, s,
    z, qrule::SingularityExtractionRule, qbuffer)

    womps = qrule.outer_quad_points

    sop = singularpart(op)
    rop = regularpart(op)

    regqrule = regularpart_quadrule(qrule)
    integrate!(rop, g, f, t, s, z, regqrule, qbuffer)

    for p in 1 : length(womps)
        x = womps[p].point
        dx = womps[p].weight

        innerintegrals!(sop, x, g, f, t, s, z, qrule, dx)
    end
end
