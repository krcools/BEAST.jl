abstract type QuasiLocalOperator <: IntegralOperator end

function allocatestorage(op::QuasiLocalOperator,
        test_functions, trial_functions,
        storage::Type{Val{:bandedstorage}})

        T = scalartype(op, test_functions, trial_functions)

        nt = Threads.nthreads()
        M = Vector{Vector{Int}}(undef, nt)
        N = Vector{Vector{Int}}(undef, nt)
        V = Vector{Vector{T}}(undef, nt)
        for i in 1:nt
            M[i] = Vector{Int}()
            N[i] = Vector{Int}()
            V[i] = Vector{T}()
        end
        # M = Int[]
        # N = Int[]
        # V = T[]

        # @show M
    
        function storeq2(v,m,n)
            tid = Threads.threadid()
            push!(M[tid],m)
            # @show length(M)
            push!(N[tid],n)
            push!(V[tid],v)
        end
    
        function freeze()
            nrows = numfunctions(test_functions)
            ncols = numfunctions(trial_functions)
            Mall = reduce(vcat, M)
            Nall = reduce(vcat, N)
            Vall = reduce(vcat, V)
            return sparse(Mall,Nall,Vall, nrows, ncols)
        end
    
        return freeze, storeq2
end


function assemblechunk!(op::QuasiLocalOperator, tfs::Space, bfs::Space, store321;
    quadstrat=defaultquadstrat(op, tfs, bfs))

    tr = assemblydata(tfs); tr == nothing && return
    br = assemblydata(bfs); br == nothing && return

    T = scalartype(op, tfs, bfs)
    tol = sqrt(eps(real(T)))

    trefs = refspace(tfs)
    brefs = refspace(bfs)

    tgeo = geometry(tfs)
    bgeo = geometry(bfs)

    num_trefs = numfunctions(trefs, domain(chart(tgeo, first(tgeo))))
    num_brefs = numfunctions(brefs, domain(chart(bgeo, first(bgeo))))

    tels, tad, tact = tr
    bels, bad, bact = br
    # @show size(tad.data)
    # @show size(bad.data)

    @assert length(tels) == size(tad.data, 3)
    @assert length(bels) == size(bad.data, 3)

    qd = quaddata(op, trefs, brefs, tels, bels, quadstrat)
    zlocal = zeros(T, num_trefs, num_brefs)
    tree = elementstree(bels)

    δ = oprange(op)
    tid = Threads.threadid()

    tid == 1 && print("dots out of 10: ")
    todo, done, pctg = length(tels), 0, 0
    for (p,(tcell,tptr)) in enumerate(zip(tels, tact))

        tc, ts = boundingbox(tcell.vertices)
        pred = (c,s) -> boxesoverlap(c,s,tc,ts + δ/2)

        for box in boxes(tree, pred)
            for q in box
                # q is the index into bels: the active basis charts
                bcell = bels[q]
                bptr = bact[q]
                @assert q <= length(bels)
                @assert q <= size(bad.data, 3)

                fill!(zlocal, 0)
                qrule = quadrule(op, trefs, brefs, p, tcell, q, bcell, qd, quadstrat)
                momintegrals!(zlocal, op,
                    tfs, tptr, tcell,
                    bfs, bptr, bcell, qrule)

                for j in 1 : length(bad[q])
                    for i in 1 : length(tad[p])
                        zij = zlocal[i,j]
                        for (n,b) in bad[q][j]
                            zb = zij*b
                            for (m,a) in tad[p][i]
                                store321(a*zb, m, n)
        end end end end end end

        done += 1
        new_pctg = round(Int, done / todo * 100)
        if new_pctg > pctg + 9
            tid == 1 && print(".")
            pctg = new_pctg
    end end
    tid == 1 && println("")
end