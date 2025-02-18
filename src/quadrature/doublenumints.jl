struct DoubleQuadRule{P,Q}
  outer_quad_points::P
  inner_quad_points::Q
end


"""
momintegrals!(biop, tshs, bshs, tcell, bcell, interactions, strat)

Function for the computation of moment integrals using simple double quadrature.
"""
function momintegrals!(biop,
    tshs, bshs, tcell, bcell, z, strat::DoubleQuadRule)

    igd = Integrand(biop, tshs, bshs, tcell, bcell)

    womps = strat.outer_quad_points
    wimps = strat.inner_quad_points
    
    for womp in womps
        tgeo = womp.point
        tvals = womp.value
        M = length(tvals)
        jx = womp.weight
        
        for wimp in wimps
            bgeo = wimp.point
            bvals = wimp.value
            N = length(bvals)
            jy = wimp.weight

            j = jx * jy

            z1 = j * igd(tgeo, bgeo, tvals, bvals)
            for n in 1:N
                for m in 1:M
                    z[m,n] += z1[m,n]
            end end
        end
    end

    return z
end


_TransposedStrat(a::DoubleQuadRule) = a