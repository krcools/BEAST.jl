abstract type SingularityExtractionStrategy end
regularpart_quadrule(qr::SingularityExtractionStrategy) = qr.regularpart_quadrule

function momintegrals!(op, g, f, t, s, z, strat::SingularityExtractionStrategy)


    womps = strat.outer_quad_points

    sop = singularpart(op)
    rop = regularpart(op)

    # compute the regular part
    rstrat = regularpart_quadrule(strat)
    momintegrals!(rop, g, f, t, s, z, rstrat)

    for p in 1 : length(womps)
        x = womps[p].point
        dx = womps[p].weight

        innerintegrals!(sop, x, g, f, t, s, z, strat, dx)
    end # next quadrature point

end
