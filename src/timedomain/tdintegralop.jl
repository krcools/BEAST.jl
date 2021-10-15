using WiltonInts84

abstract type RetardedPotential{T} <: Operator end
Base.eltype(::RetardedPotential{T}) where {T} = T
scalartype(A::RetardedPotential{T}) where {T} = T

mutable struct EmptyRP{T} <: RetardedPotential{T}
    speed_of_light::T
end
Base.eltype(::EmptyRP) = Int
defaultquadstrat(::EmptyRP, tfs, bfs) = nothing
quaddata(op::EmptyRP, xs...) = nothing
quadrule(op::EmptyRP, xs...) = nothing
momintegrals!(z, op::EmptyRP, xs...) = nothing

function allocatestorage(op::RetardedPotential, testST, basisST,
	::Type{Val{:densestorage}},
	::Type{LongDelays{:ignore}})

    tfs = spatialbasis(testST)
    bfs = spatialbasis(basisST)

    M = numfunctions(tfs)
    N = numfunctions(bfs)

    K0 = zeros(Int, M, N)
    K1 = zeros(Int, M, N)

    function store(v,m,n,k)
        K0[m,n] = (K0[m,n] == 0) ? K0[m,n] : min(K0[m,n],k)
        K1[m,n] = max(K1[m,n],k)
    end

    aux = EmptyRP(op.speed_of_light)
    print("Allocating memory for convolution operator: ")
    assemble!(aux, testST, basisST, store)
    println("\nAllocated memory for convolution operator.")

    kmax = maximum(K1);
	Z = zeros(eltype(op), M, N, kmax)
    store1(v,m,n,k) = (Z[m,n,k] += v)
    return ()->MatrixConvolution(Z), store1
end


# function allocatestorage(op::RetardedPotential, testST, basisST,
# 	::Type{Val{:bandedstorage}},
# 	::Type{LongDelays{:ignore}})

#     tfs = spatialbasis(testST)
#     bfs = spatialbasis(basisST)

#     M = numfunctions(tfs)
#     N = numfunctions(bfs)

#     K0 = fill(typemax(Int), M, N)
#     K1 = zeros(Int, M, N)

#     function store(v,m,n,k)
#         K0[m,n] = min(K0[m,n],k)
#         K1[m,n] = max(K1[m,n],k)
#     end

#     aux = EmptyRP(op.speed_of_light)
#     print("Allocating memory for convolution operator: ")
#     assemble!(aux, testST, basisST, store)
#     println("\nAllocated memory for convolution operator.")

# 	maxk1 = maximum(K1)
# 	bandwidth = maximum(K1 .- K0 .+ 1)
# 	data = zeros(eltype(op), bandwidth, M, N)
# 	Z = SparseND.Banded3D(K0, data, maxk1)
#     store1(v,m,n,k) = (Z[m,n,k] += v)
#     return ()->Z, store1
# end

struct Storage{T} end

function allocatestorage(op::RetardedPotential, testST, basisST,
	::Type{Val{:bandedstorage}},
    ::Type{LongDelays{:compress}})
    
    @info "Allocating mem for RP op compressing the static tail..."

	T = eltype(op)

    tfs = spatialbasis(testST)
    bfs = spatialbasis(basisST)

	Δt = timestep(temporalbasis(basisST))
	Nt = numfunctions(temporalbasis(basisST))

	tbf = convolve(temporalbasis(testST), temporalbasis(basisST))
	has_tail = !all(tbf.polys[end].data .== 0)

    M = numfunctions(tfs)
    N = numfunctions(bfs)

    K0 = fill(typemax(Int), M, N)
    K1 = zeros(Int, M, N)

    function store_alloc(v,m,n,k)
        K0[m,n] = min(K0[m,n],k)
        K1[m,n] = max(K1[m,n],k)
    end

    op_alloc = EmptyRP(op.speed_of_light)
	tbf_trunc = truncatetail(tbf)
	δ = timebasisdelta(Δt, Nt)
    print("Allocating memory for convolution operator: ")
    assemble!(op_alloc, tfs⊗δ, bfs⊗tbf_trunc, store_alloc)
    println("\nAllocated memory for convolution operator.")

	bandwidth = maximum(K1 .- K0 .+ 1)
	data = zeros(T, bandwidth, M, N)
    tail = zeros(T, M, N)
    # kmax = maximum(K1)
    len = has_tail ? Nt : maximum(K1)
	Z = ConvOp(data, K0, K1, tail, len)

    function store1(v,m,n,k)
		if Z.k0[m,n] ≤ k ≤ Z.k1[m,n]
			Z.data[k - Z.k0[m,n] + 1,m,n] += v
		elseif k == Z.k1[m,n]+1
			Z.tail[m,n] += v
		end
	end

    return ()->Z, store1
end

function assemble!(op::LinearCombinationOfOperators, tfs::SpaceTimeBasis, bfs::SpaceTimeBasis, store,
    threading=Threading{:multi}; quadstrat=defaultquadstrat(op, tfs, bfs))

    for (a,A) in zip(op.coeffs, op.ops)
        store1(v,m,n,k) = store(a*v,m,n,k)
        assemble!(A, tfs, bfs, store1)
    end
end

function assemble!(op::RetardedPotential, testST, trialST, store,
    threading=Threading{:multi}; quadstrat=defaultquadstrat(op, testST, trialST))

	Y, S = spatialbasis(testST), temporalbasis(testST)

	P = Threads.nthreads()
	@assert P >= 1
	splits = [round(Int,s) for s in range(0, stop=numfunctions(Y), length=P+1)]

	# @info "Starting assembly with $P threads:"
	# Threads.@threads for i in 1:P
	# 	lo, hi = splits[i]+1, splits[i+1]
	# 	lo <= hi || continue
	# 	Y_p = subset(Y, lo:hi)
	# 	store1 = (v,m,n,k) -> store(v,lo+m-1,n,k)
	# 	assemble_chunk!(op, Y_p ⊗ S, trialST, store1)
	# end

	assemble_chunk!(op, testST, trialST, store)
end

function assemble_chunk!(op::RetardedPotential, testST, trialST, store; 
    quadstrat=defaultquadstrat(op, testST, trialST))

	myid = Threads.threadid()

    testspace  = spatialbasis(testST)
    trialspace = spatialbasis(trialST)
    timebasisfunction = convolve(temporalbasis(testST), temporalbasis(trialST))

	testels, testad   = assemblydata(testspace)
	trialels, trialad = assemblydata(trialspace)

	speedoflight = op.speed_of_light
	Δt = timestep(timebasisfunction)
	ct, hs = boundingbox(vertices(geometry(trialspace)))
	diam = 2 * sqrt(3) * hs
	#kmax = ceil(Int, diam/speedoflight/timestep(timebasisfunction)) + (numintervals(timebasisfunction)-1)
	kmax = ceil(Int, (numintervals(timebasisfunction)-1) + diam/speedoflight/Δt)+1
	kmax = max(kmax, numfunctions(timebasisfunction))
    timead = temporalassemblydata(timebasisfunction, kmax=kmax)

	Δt = timestep(timebasisfunction)
	ΔR = Δt * speedoflight

	Nt = numfunctions(timebasisfunction)
	tmax = (Nt-1) * Δt

    U = refspace(testspace)
    V = refspace(trialspace)
    W = refspace(timebasisfunction)

    qd = quaddata(op, U, V, W, testels, trialels, nothing, quadstrat)

    udim = numfunctions(U)
    vdim = numfunctions(V)
    wdim = numfunctions(W)
    z = zeros(eltype(op), udim, vdim, wdim)

	@show length(testels) length(trialels)

    myid == 1 && print("dots out of 10: ")
    todo, done, pctg = length(testels), 0, 0
    for p in eachindex(testels)
        τ = testels[p]
        for q in eachindex(trialels)
            σ = trialels[q]
	        for r in rings(τ,σ,ΔR)
				r > numfunctions(timebasisfunction) && continue
	            ι = ring(r,ΔR)

	            # compute interactions between reference shape functions
	            fill!(z, 0)
	            qr = quadrule(op, U, V, W, p, τ, q, σ, r, ι, qd, quadstrat)
                momintegrals!(z, op, U, V, W, τ, σ, ι, qr)

		        # assemble in the global matrix
		        for i in 1 : udim
		            for j in 1 : vdim
		                for d in 1 : wdim

		                    v = z[i,j,d]

                            for (m,a) in testad[p,i]
                                av = a*v
                                for (n,b) in trialad[q,j]
									abv = b*av
									tad_rd = timead[r,d]
                                    for (k,c) in tad_rd
										store(c*abv, m, n, k)
									end # next κ
		                        end # next ν
		                    end # next μ
		                end # next d
		            end # next j
		        end #next i
		    end # next r
		end # next q

        done += 1
        new_pctg = round(Int, done / todo * 100)
        if myid == 1 && new_pctg > pctg + 9
            print(".")
            pctg = new_pctg
        end
    end # next p

    println("")
end


struct WiltonInts84Strat{T,V,W}
	outer_quad_points::T
    binomials::V
    workspace::W
end


function momintegrals!(z, op, g, f, T, τ, σ, ι, qr::WiltonInts84Strat)

    XW = qr.outer_quad_points
    for p in 1 : length(XW)
        x = XW[p].point
        w = XW[p].weight
		innerintegrals!(z, op, x, g, f, T, τ, σ, ι, qr, w)
    end

end
