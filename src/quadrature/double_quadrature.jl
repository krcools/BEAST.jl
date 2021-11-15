struct DoubleQuadRule{P,Q}
  outer_quad_points::P
  inner_quad_points::Q
end


"""
    regularcellcellinteractions!(biop, tshs, bshs, tcell, bcell, interactions, strat)

Function for the computation of moment integrals using simple double quadrature.
"""
function momintegrals!(biop, tshs, bshs, tcell, bcell, z, strat::DoubleQuadRule)

    # memory allocation here is a result from the type instability on strat
    # which is on purpose, i.e. the momintegrals! method is chosen based
    # on dynamic polymorphism.
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
            kernel = kernelvals(biop, tgeo, bgeo)

            for n in 1 : N
                bval = bvals[n]
                for m in 1 : M
                    tval = tvals[m]
                    igd = integrand(biop, kernel, tval, tgeo, bval, bgeo)
                    z[m,n] += j * igd
                end
            end
        end
    end

    return z
end
