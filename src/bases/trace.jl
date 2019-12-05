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

function faces(c)
    [
        simplex(c[2], c[3]),
        simplex(c[3], c[1]),
        simplex(c[1], c[2]),
    ]
end

function faces(c::CompScienceMeshes.Simplex{3,3,0,4,Float64})
    [
        simplex(c[2], c[3], c[4]),
        simplex(c[1], c[3], c[4]),
        simplex(c[1], c[2], c[4]),
        simplex(c[1], c[2], c[3])
    ]
end

function ntrace(X::Space, γ)

    x = refspace(X)
    on_target = overlap_gpredicate(γ)
    # ad = assemblydata(X)
    E, ad = assemblydata(X)

    geo = geometry(X)
    Γ = geo
    Dγ = dimension(γ)
    Σ = skeleton(Γ,Dγ)

    D = copy(transpose(connectivity(Σ, Γ, abs)))
    rows, vals = rowvals(D), nonzeros(D)

    T = scalartype(X)
    fns = [Shape{T}[] for i in 1:numfunctions(X)]

    for (p,el) in enumerate(E)

        for (q,fc) in enumerate(faces(el))

            on_target(fc) || continue
            Q = ntrace(x,el,q,fc)
            @assert norm(Q,Inf) != 0

            # find the global index in Σ of the q-th face of the p-element
            r = 0
            for k in nzrange(D,p)
                vals[k] == q && (r = rows[k]; break)
            end
            @assert r != 0
            for i in 1:size(Q,1)
                for j in 1:size(Q,2)
                    for (m,a) in ad[p,j]
                        j == q && println("bingo")
                        v = a*Q[i,j]
                        @assert a != 0
                        v == 0 && continue
                        push!(fns[m], Shape(r, i, v))
                    end
                end
            end

        end

    end

    ntrace(X, Σ, fns)
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
            cell = Σ.faces[e]
            fc = chart(Σ, cell)
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

"""
currently not working!
"""
function ttrace(X::Space, γ)

    x = refspace(X)
    on_target = overlap_gpredicate(γ)
    # ad = assemblydata(X)
    E, ad = assemblydata(X)

    geo = geometry(X)
    Γ = geo
    Dγ = dimension(γ)
    Σ = skeleton(Γ,Dγ)

    D = copy(transpose(connectivity(Σ, Γ, abs)))
    rows, vals = rowvals(D), nonzeros(D)

    T = scalartype(X)
    fns = [Shape{T}[] for i in 1:numfunctions(X)]

    for (p,el) in enumerate(E)

        for (q,fc) in enumerate(faces(el))

            on_target(fc) || continue
            Q = ttrace(x,el,q,fc)
            print(Q, "\n")
            # find the global index in Σ of the q-th face of the p-element
            r = 0
            for k in nzrange(D,p)
                vals[k] == q && (r = rows[k]; break)
            end
            @assert r != 0
            for i in 1:size(Q,1)
                for j in 1:size(Q,2)
                    for (m,a) in ad[p,j]
                        if j ≈ q
                            print("bingo")
                        end
                        v = a*Q[i,j]
                        v == 0 && continue
                        push!(fns[m], Shape(r, i, v))
                    end
                end
            end

        end

    end

    ttrace(X, Σ, fns)
end
