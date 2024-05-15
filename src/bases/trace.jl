using SparseArrays

"""
    ntrace(X::Space, γ::Mesh)

Compute the normal trace of basis X on mesh γ. γ is assumed to be part of the
boundary of geometry(X).
"""
function ntrace(X::DirectProductSpace{T}, γ) where T
    x = Space{T}[ntrace(s, γ) for s in X.factors]
    return DirectProductSpace(x)
end


function strace(X::DirectProductSpace{T}, γ) where T
    x = Space{T}[strace(s,γ) for s in X.factors]
    return DirectProductSpace(x)
end



function ntrace(X::Space, γ)

    # on_target = overlap_gpredicate(γ)
    # ad = assemblydata(X)
    x = refspace(X)
    E, ad, P = assemblydata(X)
    igeo = geometry(X)
    @assert dimension(γ) == dimension(igeo)-1
    # Γ = geo
    # Dγ = dimension(γ)
    # Σ = skeleton(Γ,Dγ)

    ogeo = boundary(igeo)
    on_target = overlap_gpredicate(γ)
    ogeo = submesh(ogeo) do m,f
        ch = chart(m,f)
        on_target(ch)
    end

    # D = copy(transpose(connectivity(ogeo, igeo, abs)))
    D = connectivity(igeo, ogeo, abs)
    rows, vals = rowvals(D), nonzeros(D)

    T = scalartype(X)
    S = Shape{T}
    fns = [Vector{S}() for i in 1:numfunctions(X)]

    for (p,el) in enumerate(E)

        for (q,fc) in enumerate(faces(el))
            on_target(fc) || continue

            # print(Q)
            # @assert norm(Q,Inf) != 0
            
            r = 0
            for k in nzrange(D,P[p])
                vals[k] == q && (r = rows[k]; break)
            end
            @assert r != 0
            
            fc1 = chart(ogeo, r)
            Q = ntrace(x, el, q, fc1)

            for i in 1:size(Q,1)
                for j in 1:size(Q,2)
                    for (m,a) in ad[p,j]
                        # j == q && println("bingo",j,q)
                        v = a*Q[i,j]
                        isapprox(v,0,atol=sqrt(eps(T))) && continue
                        push!(fns[m], Shape(r, i, v))
                    end
                end
            end

        end

    end

    ntrace(X, ogeo, fns)
end


strace(X::Space, γ::Mesh) = strace(X, γ, Val{dimension(γ)+1})

function strace(X::Space, γ, dim1::Type{Val{2}})

    x = refspace(X)
    on_target = overlap_gpredicate(γ)
    E, ad = assemblydata(X)

    geo = geometry(X)
    Γ = geo
    Σ = skeleton(Γ,1)

    # Allows for quick retrieval of the faces in Σ in a given cell of Γ.
    D = copy(transpose(connectivity(Σ, Γ, abs)))
    rows, vals = rowvals(D), nonzeros(D)

    T = scalartype(X)
    @assert numfunctions(X) != 0
    fns = [Shape{T}[] for i in 1:numfunctions(X)]

    for (p,el) in enumerate(E)

        for q in 1:dimension(Γ)+1

            # find the global index in Σ of the q-th face of the p-element
            e = 0
            for k in nzrange(D,p)
                if vals[k] == q
                    e = rows[k]
                    break
                end
            end
            @assert e != 0

            # make sure we use fc as oriented in Σ
            cell = CompScienceMeshes.indices(Σ,e) #Σ.faces[e]
            fc = chart(Σ, e)
            on_target(fc) || continue

            Q = strace(x,el,q,fc)
            for i in 1:size(Q,1)
                for j in 1:size(Q,2)
                    for (m,a) in ad[p,j]
                        v = a*Q[i,j]
                        v == 0 && continue
                        push!(fns[m], Shape(e, i, v))
                    end
                end
            end

        end

    end

    strace(X, Σ, fns)
end

function strace(X::Space, γ, dim1::Type{Val{3}})

    x = refspace(X)
    E, ad, P = assemblydata(X)
    igeo = geometry(X)

    @assert dimension(γ) == dimension(igeo)-1

    ogeo = boundary(igeo)
    on_target = overlap_gpredicate(γ)
    ogeo = submesh(ogeo) do m,f
        ch = chart(m,f)
        on_target(ch)
    end

    D = connectivity(igeo, ogeo, abs)
    rows, vals = rowvals(D), nonzeros(D)

    T = scalartype(X)
    S = Shape{T}
    fns = [Vector{S}() for i in 1:numfunctions(X)]

    for (p,el) in enumerate(E)

        for (q,fc) in enumerate(faces(el))
            on_target(fc) || continue

            r = 0
            for k in nzrange(D,P[p])
                vals[k] == q && (r = rows[k]; break)
            end
            @assert r != 0

            fc1 = chart(ogeo, r)
            Q = strace(x, el, q, fc1)

            for i in 1:size(Q,1)
                for j in 1:size(Q,2)
                    for (m,a) in ad[p,j]
                        v = a*Q[i,j]
                        isapprox(v,0,atol=sqrt(eps(T))) && continue
                        push!(fns[m], Shape(r, i, v))
                    end
                end
            end

        end

    end

    strace(X, ogeo, fns)

end

"""
currently not working!
"""
# function ttrace(X::Space, γ)
#
#     x = refspace(X)
#     on_target = overlap_gpredicate(γ)
#     # ad = assemblydata(X)
#     E, ad = assemblydata(X)
#
#     geo = geometry(X)
#     Γ = geo
#     Dγ = dimension(γ)
#     Σ = skeleton(Γ,Dγ)
#
#     D = copy(transpose(connectivity(Σ, Γ, abs)))
#     rows, vals = rowvals(D), nonzeros(D)
#
#     T = scalartype(X)
#     fns = [Shape{T}[] for i in 1:numfunctions(X)]
#
#     for (p,el) in enumerate(E)
#
#         for (q,fc) in enumerate(faces(el))
#
#             on_target(fc) || (println("skip"); continue)
#             # print("\n")
#             Q = ttrace(x,el,q,fc)
#             # print(Q,q)
#             # find the global index in Σ of the q-th face of the p-element
#             r = 0
#             for k in nzrange(D,p)
#                 vals[k] == q && (r = rows[k]; break)
#             end
#             @assert r != 0
#             for i in 1:size(Q,1)
#                 for j in 1:size(Q,2)
#                     for (m,a) in ad[p,j]
#                         # if j == q
#                         #     # print("bingo",a,Q[i,j])
#                         # end
#                         v = a*Q[i,j]
#                         v == 0 && continue
#                         push!(fns[m], Shape(r, i, v))
#                     end
#                 end
#             end
#
#         end
#
#     end
#
#     ttrace(X, Σ, fns)
# end


function ttrace(X::Space, γ)

    x = refspace(X)
    E, ad = assemblydata(X)
    igeo = geometry(X)
    @assert dimension(γ) == dimension(igeo)-1

    ogeo = boundary(igeo)
    on_target = overlap_gpredicate(γ)
    ogeo = submesh(ogeo) do m,face
        ch = chart(m,face)
        on_target(ch)
    end

    D = connectivity(igeo, ogeo, abs)
    rows, vals = rowvals(D), nonzeros(D)

    T = scalartype(X)
    S = Shape{T}
    fns = [Vector{S}() for i in 1:numfunctions(X)]

    for (p,el) in enumerate(E)

        for (q,fc) in enumerate(faces(el))
            on_target(fc) || continue

            r = 0
            for k in nzrange(D,p)
                vals[k] == q && (r = rows[k]; break)
            end
            @assert r != 0

            fc1 = chart(ogeo, r)
            @assert cartesian(center(fc)) ≈ cartesian(center(fc1))

            Q = ttrace(x, el, q, fc1)
            for i in 1:size(Q,1)
                for j in 1:size(Q,2)
                    for (m,a) in ad[p,j]
                        v = a*Q[i,j]
                        isapprox(v,0,atol=sqrt(eps(T))) && continue
                        push!(fns[m], Shape(r, i, v))
                    end
                end
            end

        end
    end

    ttrace(X, ogeo, fns)
end
