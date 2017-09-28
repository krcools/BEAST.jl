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
    zlocal = zeros(scalartype(biop, tfs, bfs), num_tshapes, num_bshapes)

    assemblechunk_body!(biop,
        tshapes, test_elements, tad,
        bshapes, bsis_elements, bad,
        qd, zlocal, store)
end

function assemblerow!(biop::IntegralOperator, test_functions::Space, trial_functions::Space, store)

    test_elements = elements(geometry(test_functions))
    trial_elements, trial_assembly_data = assemblydata(trial_functions)

    test_shapes = refspace(test_functions); num_test_shapes = numfunctions(test_shapes)
    trial_shapes = refspace(trial_functions); num_trial_shapes = numfunctions(trial_shapes)

    qd = quaddata(biop, test_shapes, trial_shapes, test_elements, trial_elements)
    zlocal = zeros(scalartype(biop, test_functions, trial_functions), num_test_shapes, num_trial_shapes)

    @assert size(zlocal) == (3,3)

    @assert length(trial_elements) == numcells(geometry(trial_functions))
    @assert numfunctions(test_functions) == 1
    test_function = test_functions.fns[1]
    for shape in test_function
        p = shape.cellid
        i = shape.refid
        a = shape.coeff
        #@show shape
        tcell = test_elements[p]
        for (q,bcell) in enumerate(trial_elements)

            #@show q

            fill!(zlocal, 0)
            strat = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, qd)
            momintegrals!(biop, test_shapes, trial_shapes, tcell, bcell, zlocal, strat)

            #@assert size(zlocal,2) == 3
            for j in 1:size(zlocal,2)
                for (n,b) in trial_assembly_data[q,j]
                    #n == 3 && println("bingo")
                    #@show n
                    store(a*zlocal[i,j]*b, 1, n)
                    #store(1, 1, n)
end end end end end


# function assemblerow_body!(biop::IntegralOperator,
#     tshapes, test_function,
#     bshapes, trial_elements, trial_assembly_data,
#     qd, zlocal, store)
#
#     for shape in test_function
#         p = shape.cellid
#         tcell = test_elements[p]
#         for (q,bcell) in enumerate(trial_elements)
#
#             fill!(zlocal, 0)
#             strat = quadrule(biop, tshapes, bshapes, p, tcell, q, bcell, qd)
#             momintegrals!(biop, test_shapes, trial_shapes, tcell, bcell, zlocal, strat)
#
#             i = shape.refid
#             a = shape.coeff
#             for j in size(zlocal,2)
#                 for (n,b) in trial_assembly_data[q,j]
#                     store(a*zlocal[i,j]*b, m, n)
# end end end end end



function assemblechunk_body!(biop,
        test_shapes, test_elements, test_assembly_data,
        trial_shapes, trial_elements, trial_assembly_data,
        qd, zlocal, store)

    for (p,tcell) in enumerate(test_elements), (q,bcell) in enumerate(trial_elements)

        fill!(zlocal, 0)
        strat = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, qd)
        momintegrals!(biop, test_shapes, trial_shapes, tcell, bcell, zlocal, strat)

        for j in 1 : size(zlocal,2), i in 1 : size(zlocal,1)
            for (n,b) in trial_assembly_data[q,j], (m,a) in test_assembly_data[p,i]
                store(a*zlocal[i,j]*b, m, n)
end end end end






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
