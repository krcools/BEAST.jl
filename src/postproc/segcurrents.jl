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
    Hi
"""
function grideval(points, coeffs, basis)

    # charts: active charts
    # ad: assembly data (active_cell_idx, local_shape_idx) -> [dof1, dfo2, ...]
    # ag: active_cell_idx -> global_cell_idx
    charts, ad, ag = assemblydata(basis)
    refs = refspace(basis)

    V = valuetype(refs, eltype(charts))
    T = promote_type(eltype(coeffs), eltype(V))
    P = similar_type(V, T)

    values = zeros(P, length(points))

    chart_tree = BEAST.octree(charts)
    for (j,point) in enumerate(points)
        i = CompScienceMeshes.findchart(charts, chart_tree, point)
        if i != nothing
            chart = charts[i]
            u = carttobary(chart, point)
            vals = refs(neighborhood(chart,u))
            for r in 1 : numfunctions(refs)
                for (m,w) in ad[i, r]
                    values[j] += w * coeffs[m] * vals[r][1]
                end
            end
        end
    end
    return values
end
