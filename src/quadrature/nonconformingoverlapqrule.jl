struct NonConformingOverlapQRule{S}
    conforming_qstrat::S
end


function momintegrals!(op,
    test_local_space, basis_local_space,
    test_chart::CompScienceMeshes.Simplex, basis_chart::CompScienceMeshes.Simplex,
    out, qrule::NonConformingOverlapQRule)

    test_charts, tclps = CompScienceMeshes.intersection_keep_clippings(test_chart, basis_chart)
    _, bclps = CompScienceMeshes.intersection_keep_clippings(basis_chart, test_chart)
    bsis_charts = copy(test_charts)

    for tclp in tclps append!(test_charts, tclp) end
    for bclp in bclps append!(bsis_charts, bclp) end

    T = coordtype(test_chart)
    h = max(volume(test_chart), volume(test_chart))
    test_charts = [ch for ch in test_charts if volume(ch) .> 1e6 * eps(T) * h]
    bsis_charts = [ch for ch in bsis_charts if volume(ch) .> 1e6 * eps(T) * h]

    test_charts = [ch for ch in test_charts if volume(ch) .> 1e6 * eps(T)]
    bsis_charts = [ch for ch in bsis_charts if volume(ch) .> 1e6 * eps(T)]

    isempty(test_charts) && return
    isempty(bsis_charts) && return

    # @assert volume(test_chart) ≈ sum(volume.(test_charts))
    # if volume(basis_chart) ≈ sum(volume.(bsis_charts)) else
    #     @show volume(basis_chart)
    #     @show sum(volume.(bsis_charts))
    #     error()
    # end
    # test_local_space = refspace(test_functions)
    # basis_local_space = refspace(basis_functions)

    qstrat = CommonFaceOverlappingEdgeQStrat(qrule.conforming_qstrat)
    qdata = quaddata(op, test_local_space, basis_local_space,
        test_charts, bsis_charts, qstrat)

    for (p,tchart) in enumerate(test_charts)
        for (q,bchart) in enumerate(bsis_charts)
            qrule = quadrule(op, test_local_space, basis_local_space,
                p, tchart, q, bchart, qdata, qstrat)
            # @show qrule

            P = restrict(test_local_space, test_chart, tchart)
            Q = restrict(basis_local_space, basis_chart, bchart)
            zlocal = zero(out)

            # momintegrals!(zlocal, op,
            #     test_local_space, nothing, tchart,
            #     basis_local_space, nothing, bchart, qrule)
            momintegrals!(op, test_local_space, basis_local_space,
                tchart, bchart, zlocal, qrule)

            for i in axes(P,1)
                for j in axes(Q,1)
                    for k in axes(P,2)
                        for l in axes(Q,2)
                            out[i,j] += P[i,k] * zlocal[k,l] * Q[j,l]
            # out .+= P * zlocal * Q'    
end end end end end end end