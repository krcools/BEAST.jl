struct NonConformingTouchQRule{S}
    conforming_qstrat::S
    test_overlapping_edge_index::Int
    bsis_overlapping_edge_index::Int
end

function momintegrals!(op,
    test_locspace, bsis_locspace,
    τ::CompScienceMeshes.Simplex, σ::CompScienceMeshes.Simplex,
    out, qrule::NonConformingTouchQRule)

    num_tshapes = numfunctions(test_locspace, domain(τ))
    num_bshapes = numfunctions(bsis_locspace, domain(σ))

    T = coordtype(τ)
    P = eltype(τ.vertices)

    i = qrule.test_overlapping_edge_index
    j = qrule.bsis_overlapping_edge_index

    @assert volume(τ) > eps(T) * 1e3
    @assert volume(σ) > eps(T) * 1e3
    τs, σs = _conforming_refinement_touching_triangles(τ,σ,i,j)

    isempty(τs) && return
    isempty(σs) && return

    @assert all(volume.(τs) .> 1e3 * eps(T) * (volume(τ)))
    @assert all(volume.(σs) .> 1e3 * eps(T) * (volume(σ)))
 
    test_overlaps = map(τs) do tchart
        simplex(
            carttobary(τ, tchart.vertices[1]),
            carttobary(τ, tchart.vertices[2]),
            carttobary(τ, tchart.vertices[3]))
    end

    trial_overlaps = map(σs) do bchart
        simplex(
            carttobary(σ, bchart.vertices[1]),
            carttobary(σ, bchart.vertices[2]),
            carttobary(σ, bchart.vertices[3]))
    end

    qstrat = qrule.conforming_qstrat
    qdata = quaddata(op, test_locspace, bsis_locspace, τs, σs, qstrat)

    @assert !any(volume.(τs) .< 1e-13)
    @assert !any(volume.(σs) .< 1e-13)

    zlocal = zero(out)
    P = zeros(T, num_tshapes, num_tshapes)
    Q = zeros(T, num_bshapes, num_bshapes)
    for (p,tchart) in enumerate(τs)
        restrict!(P, test_locspace, τ, tchart, test_overlaps[p])
        for (q,bchart) in enumerate(σs)
            restrict!(Q, bsis_locspace, σ, bchart, trial_overlaps[q])
            
            qrule = quadrule(op, test_locspace, bsis_locspace,
                p, tchart, q, bchart, qdata, qstrat)

            fill!(zlocal, 0)
            momintegrals!(op, test_locspace, bsis_locspace,
                tchart, bchart, zlocal, qrule)

            for i in axes(P,1)
                for j in axes(Q,1)
                    for k in axes(P,2)
                        for l in axes(Q,2)
                            out[i,j] += P[i,k] * zlocal[k,l] * Q[j,l]
end end end end end end end


function _conforming_refinement_touching_triangles_bak(τ,σ,i,j)
    λ = faces(τ)[i]
    μ = faces(σ)[j]
    ρ = CompScienceMeshes.intersection(λ,μ)[1]

    τ_verts = [τ[mod1(i+2,3)], τ[mod1(i,3)], τ[mod1(i+1,3)]]
    σ_verts = [σ[mod1(j+2,3)], σ[mod1(j,3)], σ[mod1(j+1,3)]]

    T = coordtype(τ)
    P = eltype(τ.vertices)

    U = T[]
    V = P[]
    for v in ρ.vertices
        if CompScienceMeshes.isinside(λ,v)
            push!(V,v)
            push!(U, carttobary(λ,v)[1])
    end end
    if length(U) == 2
        if U[1] < U[2]
            temp = V[1]
            V[1] = V[2]
            V[2] = temp
    end end
    append!(τ_verts, V)

    U = T[]
    V = P[]
    for v in ρ.vertices
        if CompScienceMeshes.isinside(μ,v)
            push!(V,v)
            push!(U, carttobary(μ,v)[1])
    end end
    if length(U) == 2
        if U[1] < U[2]
            temp = V[1]
            V[1] = V[2]
            V[2] = temp
    end end
    append!(σ_verts, V)

    τ_verts = push!(τ_verts[2:end], τ_verts[1])
    σ_verts = push!(σ_verts[2:end], σ_verts[1])

    # @show τ_verts

    τ_charts = [ simplex(τ_verts[1], τ_verts[i], τ_verts[i+1]) for i in 2:length(τ_verts)-1 ]
    σ_charts = [ simplex(σ_verts[1], σ_verts[i], σ_verts[i+1]) for i in 2:length(σ_verts)-1 ]

    signs = Int.(sign.(dot.(normal.(τ_charts),Ref(normal(τ)))))
    τ_charts = flip_normal.(τ_charts,signs)
    signs = Int.(sign.(dot.(normal.(σ_charts),Ref(normal(σ)))))
    σ_charts = flip_normal.(σ_charts,signs)

    h = sqrt(volume(τ))
    τ_charts = τ_charts[volume.(τ_charts) .> 1e3 * eps(T) * h]
    σ_charts = τ_charts[volume.(σ_charts) .> 1e3 * eps(T) * h]

    return τ_charts, σ_charts
end


function _conforming_refinement_touching_triangles(τ,σ,i,j)
    λ = faces(τ)[i]
    μ = faces(σ)[j]
    ρ = CompScienceMeshes.intersection(λ,μ)[1]

    τ_verts = Array(τ.vertices)
    σ_verts = Array(σ.vertices)

    T = coordtype(τ)
    P = eltype(τ.vertices)

    U = T[]
    V = P[]
    for v in ρ.vertices
        push!(V,v)
        push!(U, carttobary(λ,v)[1])
    end
    if U[1] < U[2]
        temp = V[1]
        V[1] = V[2]
        V[2] = temp
    end
    p = mod1(i+2,3)
    new_i = p <= i ? i+2 : i
    insert!(τ_verts, p, V[2])
    insert!(τ_verts, p, V[1])

    U = T[]
    V = P[]
    for v in ρ.vertices
        push!(V,v)
        push!(U, carttobary(μ,v)[1])
    end
    if U[1] < U[2]
        temp = V[1]
        V[1] = V[2]
        V[2] = temp
    end
    p = mod1(j+2,3)
    new_j = p <= j ? j+2 : j
    insert!(σ_verts, p, V[2])
    insert!(σ_verts, p, V[1])

    τ_charts = [ simplex(τ_verts[mod1(new_i,5)], τ_verts[mod1(new_i+s,5)], τ_verts[mod1(new_i+s+1,5)]) for s in 1:3 ]
    σ_charts = [ simplex(σ_verts[mod1(new_j,5)], σ_verts[mod1(new_j+s,5)], σ_verts[mod1(new_j+s+1,5)]) for s in 1:3 ]

    τ_charts = [ch for ch in τ_charts if volume(ch) .> 1e6 * eps(T) * volume(τ)]
    σ_charts = [ch for ch in σ_charts if volume(ch) .> 1e6 * eps(T) * volume(σ)]
    τ_charts = [ch for ch in τ_charts if volume(ch) .> 1e6 * eps(T)]
    σ_charts = [ch for ch in σ_charts if volume(ch) .> 1e6 * eps(T)]

    # signs = Int.(sign.(dot.(normal.(τ_charts),Ref(normal(τ)))))
    # @assert all(signs .== 1)
    # τ_charts = flip_normal.(τ_charts,signs)
    # signs = Int.(sign.(dot.(normal.(σ_charts),Ref(normal(σ)))))
    # @assert all(signs .== 1)
    # σ_charts = flip_normal.(σ_charts,signs)

    return τ_charts, σ_charts
end


function _num_common_vertices(τ, σ)
    hits = 0
    T = coordtype(σ)
    tol = eps(T) * 10^3
    for v in τ.vertices
        for w in σ.vertices
            if norm(v - w) < tol
                hits += 1
                break
    end end end
    return hits
end


function _test_conformity(τ, σ)
    if CompScienceMeshes.overlap(τ, σ)
        # if _num_common_vertices(τ, σ) != 3
        #     @infiltrate
        # end
        return _num_common_vertices(τ, σ) == 3
    end

    for eτ in faces(τ)
        for eσ in faces(σ)
            if CompScienceMeshes.overlap(eτ, eσ)
                # if _num_common_vertices(τ, σ) != 2
                #     @infiltrate
                # end
                return _num_common_vertices(τ, σ) == 2
            end
    end end

    for u in τ.vertices
        for v in σ.vertices
            su = simplex([u], Val{0})
            sv = simplex([v], Val{0})
            if CompScienceMeshes.overlap(su, sv)
                # if _num_common_vertices(τ, σ) != 1
                #     @infiltrate
                # end
                return _num_common_vertices(τ, σ) == 1
    end end end

    @assert _num_common_vertices(τ, σ) == 0
    return _num_common_vertices(τ, σ) == 0
end