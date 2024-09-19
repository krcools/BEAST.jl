struct NonConformingOverlapQRule{S}
    conforming_qstrat::S
end
# function tangent_rank(p::CompScienceMeshes.Simplex{U,D}) where {U,D}
#     G = [dot(p.tangents[i], p.tangents[j]) for i in 1:D, j in 1:D]
#     return rank(G) == D 
# end

function momintegrals!(op,
    test_local_space, basis_local_space,
    test_chart::CompScienceMeshes.Simplex, basis_chart::CompScienceMeshes.Simplex,
    out, qrule::NonConformingOverlapQRule)

    num_tshapes = numfunctions(test_local_space, domain(test_chart))
    num_bshapes = numfunctions(basis_local_space, domain(basis_chart))

    test_charts, tclps = CompScienceMeshes.intersection_keep_clippings(test_chart, basis_chart)
    _,           bclps = CompScienceMeshes.intersection_keep_clippings(basis_chart, test_chart)
    bsis_charts = copy(test_charts)

    for tclp in tclps append!(test_charts, tclp) end
    for bclp in bclps append!(bsis_charts, bclp) end

    T = coordtype(test_chart)
    h = max(volume(test_chart), volume(test_chart))
    test_charts = [ch for ch in test_charts if volume(ch) .> 1e6 * eps(T) * h]
    bsis_charts = [ch for ch in bsis_charts if volume(ch) .> 1e6 * eps(T) * h]

    test_charts = [ch for ch in test_charts if volume(ch) .> 1e6 * eps(T)]
    bsis_charts = [ch for ch in bsis_charts if volume(ch) .> 1e6 * eps(T)]

    # test_charts = [ch for ch in test_charts if tangent_rank(ch)]
    # bsis_charts = [ch for ch in bsis_charts if tangent_rank(ch)]

    isempty(test_charts) && return
    isempty(bsis_charts) && return

    test_overlaps = map(test_charts) do tchart
        simplex(
            carttobary(test_chart, tchart.vertices[1]),
            carttobary(test_chart, tchart.vertices[2]),
            carttobary(test_chart, tchart.vertices[3]))
    end

    trial_overlaps = map(bsis_charts) do bchart
        simplex(
            carttobary(basis_chart, bchart.vertices[1]),
            carttobary(basis_chart, bchart.vertices[2]),
            carttobary(basis_chart, bchart.vertices[3]))
    end

    qstrat = CommonFaceOverlappingEdgeQStrat(qrule.conforming_qstrat)
    qdata = quaddata(op, test_local_space, basis_local_space,
        test_charts, bsis_charts, qstrat)

    zlocal = zero(out)
    P = zeros(T, num_tshapes, num_tshapes)
    Q = zeros(T, num_bshapes, num_bshapes)
    for (p,tchart) in enumerate(test_charts)
        restrict!(P, test_local_space, test_chart, tchart, test_overlaps[p])
        for (q,bchart) in enumerate(bsis_charts)
            restrict!(Q, basis_local_space, basis_chart, bchart, trial_overlaps[q])

            qrule = quadrule(op, test_local_space, basis_local_space,
                p, tchart, q, bchart, qdata, qstrat)

            fill!(zlocal, 0)
            momintegrals!(op, test_local_space, basis_local_space,
                tchart, bchart, zlocal, qrule)

            for i in axes(P,1)
                for j in axes(Q,1)
                    for k in axes(P,2)
                        for l in axes(Q,2)
                            out[i,j] += P[i,k] * zlocal[k,l] * Q[j,l]
end end end end end end end