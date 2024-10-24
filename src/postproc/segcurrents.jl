function octree(charts::Vector{S} where {S <: CompScienceMeshes.Simplex})

    ncells = length(charts)

    T = coordtype(charts[1])
    P = eltype(charts[1].vertices)

    points = zeros(P, ncells)
    radii = zeros(T, ncells)

    for (i,ch) in enumerate(charts)
        points[i] = cartesian(center(ch))
        radii[i] = maximum(norm(v-points[i]) for v in ch.vertices)
    end

    return Octree(points, radii)
end

"""
    grideval(points, coeffs, basis; type=nothing)
"""
function grideval(points, coeffs, basis; type=nothing)

    # charts: active charts
    # ad: assembly data (active_cell_idx, local_shape_idx) -> [dof1, dfo2, ...]
    # ag: active_cell_idx -> global_cell_idx
    charts, ad, ag = assemblydata(basis)
    refs = refspace(basis)

    V = valuetype(refs, eltype(charts))
    T = promote_type(eltype(coeffs), eltype(V))
    if !(V <: SVector)
        P = T
    else
        P = similar_type(V, T)
    end

    type !== nothing && (P = type)

    values = zeros(P, size(points))

    chart_tree = BEAST.octree(charts)
    for (j,point) in enumerate(points)
        i = CompScienceMeshes.findchart(charts, chart_tree, point)
        if i !== nothing
            # @show i
            chart = charts[i]
            u = carttobary(chart, point)
            vals = refs(neighborhood(chart,u))
            for r in 1 : numfunctions(refs, domain(chart))
                for (m,w) in ad[i, r]
                    values[j] += w * coeffs[m] * vals[r][1]
                end
            end
            continue
        end
    end
    return values
end
