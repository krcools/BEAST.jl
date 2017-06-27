abstract type IntegralOperator <: Operator end

export quadrule, elements

"""
    quadrule(operator,test_refspace,trial_refspace,p,test_element,q_trial_element, qd)

Returns an object that contains all the dynamic (runtime) information that
defines the integration strategy that will be used by `momintegrals!` to compute
the interactions between the local test/trial functions defined on the specified
geometric elements. The indices `p` and `q` refer to the position of the test
and trial elements as encountered during iteration over the output of
`geometry`.

The last argument `qd` provides access to all precomputed data required for
quadrature. For example it might be desirable to precompute all the quadrature
points for all possible numerical quadrature schemes that can potentially be
required during matrix assembly. This makes sense, since the number of point is
order N (where N is the number of faces) but these points will appear in N^2
computations. Precomputation requires some extra memory but can save a lot on
computation time.
"""
function quadrule end


"""
  elements(geo)

Create an iterable collection of the elements stored in `geo`. The order in which
this collection produces the elements determines the index used for lookup in the
data structures returned by `assemblydata` and `quaddata`.
"""
#elements(geo) = [simplex(vertices(geo,cl)) for cl in cells(geo)]
elements(geo) = [chart(geo,cl) for cl in cells(geo)]

elements(sp::Space) = elements(geometry(sp))

"""
    assemblechunk!(biop::IntegralOperator, tfs, bfs, store)

Computes the matrix of operator biop wrt the finite element spaces tfs and bfs
"""
function assemblechunk!(biop::IntegralOperator, tfs::Space, bfs::Space, store)

    test_elements, tad = assemblydata(tfs)
    bsis_elements, bad = assemblydata(bfs)

    tshapes = refspace(tfs); num_tshapes = numfunctions(tshapes)
    bshapes = refspace(bfs); num_bshapes = numfunctions(bshapes)

    qd = quaddata(biop, tshapes, bshapes, test_elements, bsis_elements)
    T = promote_type(scalartype(biop), scalartype(tfs), scalartype(bfs))
    zlocal = zeros(T, num_tshapes, num_bshapes)

    print("dots out of 10: ")
    todo, done, pctg = length(test_elements), 0, 0
    for p in eachindex(test_elements)
        tcell = test_elements[p]
        for q in eachindex(bsis_elements)
            bcell = bsis_elements[q]

            fill!(zlocal, 0)
            strat = quadrule(biop, tshapes, bshapes, p, tcell, q, bcell, qd)
            momintegrals!(biop, tshapes, bshapes, tcell, bcell, zlocal, strat)

            for j in 1 : num_bshapes, i in 1 : num_tshapes
                z = zlocal[i,j]
                for (n,b) in bad[q,j], (m,a) in tad[p,i]
                    store(a*z*b, m, n)
        end end end

        done += 1
        new_pctg = round(Int, done / todo * 100)
        (new_pctg > pctg + 9) && (print("."); pctg = new_pctg)
    end
    print(" done. ")
end








immutable DoubleQuadStrategy{P,Q}
  outer_quad_points::P
  inner_quad_points::Q
end


"""
    regularcellcellinteractions!(biop, tshs, bshs, tcell, bcell, interactions, strat)

Function for the computation of moment integrals using simple double quadrature.
"""
function momintegrals!(biop, tshs, bshs, tcell, bcell, z, strat::DoubleQuadStrategy)

    # memory allocation here is a result from the type instability on strat
    # which is on purpose, i.e. the momintegrals! method is chosen based
    # on dynamic polymorphism.
    womps = strat.outer_quad_points
    wimps = strat.inner_quad_points

    M, N = size(z)

    for womp in womps
        tgeo = womp.point
        tvals = womp.value
        jx = womp.weight

        for wimp in wimps
            bgeo = wimp.point
            bvals = wimp.value
            jy = wimp.weight

            j = jx * jy
            kernel = kernelvals(biop, tgeo, bgeo)

            for m in 1 : M
                tval = tvals[m]
                for n in 1 : N
                    bval = bvals[n]

                    igd = integrand(biop, kernel, tval, tgeo, bval, bgeo)
                    z[m,n] += j * igd
                end
            end
        end
    end

    return z
end


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


type QuadData{WPV1,WPV2}
  tpoints::Matrix{Vector{WPV1}}
  bpoints::Matrix{Vector{WPV2}}
end
