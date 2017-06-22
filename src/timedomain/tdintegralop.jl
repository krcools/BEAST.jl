using WiltonInts84

abstract type RetardedPotential{T} <: Operator end
Base.eltype{T}(::RetardedPotential{T}) = T
scalartype{T}(A::RetardedPotential{T}) = T

type EmptyRP{T} <: RetardedPotential{T}
    speed_of_light::T
end
Base.eltype(::EmptyRP) = Int
quaddata(op::EmptyRP, xs...) = nothing
quadrule(op::EmptyRP, xs...) = nothing
momintegrals!(z, op::EmptyRP, xs...) = nothing

function allocatestorage(op::RetardedPotential, testST, basisST)

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
    println("Allocating memory for convolution operator....")
    assemble!(aux, testST, basisST, store)
    println("Allocated memory for convolution operator.")
    data = zeros(eltype(op), M, N, maximum(K1-K0+1))

    Z = SparseND.Banded3D(K0,K1,data)
    store1(v,m,n,k) = (Z[m,n,k] += v)
    return Z, store1
end


function assemble!(op::LinearCombinationOfOperators, tfs::SpaceTimeBasis, bfs::SpaceTimeBasis, store)
    for (a,A) in zip(op.coeffs, op.ops)
        store1(v,m,n,k) = store(a*v,m,n,k)
        assemble!(A, tfs, bfs, store1)
    end
end


function assemble!(op::RetardedPotential, testST, trialST, store)

    testspace  = spatialbasis(testST)
    trialspace = spatialbasis(trialST)
    timebasisfunction = convolve(temporalbasis(testST), temporalbasis(trialST))

	testels, testad   = assemblydata(testspace)
	trialels, trialad = assemblydata(trialspace)

    timeels, timead = assemblydata(timebasisfunction)

	Δt = timestep(timebasisfunction)
	ΔR = Δt * op.speed_of_light
	Nt = numfunctions(timebasisfunction)
	tmax = (Nt-1) * Δt

    U = refspace(testspace)
    V = refspace(trialspace)
    W = refspace(timebasisfunction)

    qd = quaddata(op, U, V, W, testels, trialels, timeels)

    udim = numfunctions(U)
    vdim = numfunctions(V)
    wdim = numfunctions(W)
    z = zeros(eltype(op), udim, vdim, wdim)

    print("dots out of 10: ")
    todo, done, pctg = length(testels), 0, 0
    for p in eachindex(testels)
        τ = testels[p]
        for q in eachindex(trialels)
            σ = trialels[q]
	        for r in rings(τ,σ,ΔR)
				r > numfunctions(timebasisfunction) && continue
	            # construct the radial interval [(r-1)ΔR, rΔR]
	            ι = ring(r,ΔR)

	            # compute interactions between reference shape functions
	            fill!(z, 0)
	            qr = quadrule(op, U, V, W, p, τ, q, σ, r, ι, qd)
				momintegrals!(z, op, U, V, W, τ, σ, tmax, ι, qr)

		        # assemble in the global matrix
		        for i in 1 : udim
		            for j in 1 : vdim
		                for d in 1 : wdim

		                    v = z[i,j,d]

                            for (m,a) in testad[p,i]
                                for (n,b) in trialad[q,j]
                                    for (k,c) in timead[Nt-r,d]
										store(a*b*c*v, m, n, Nt-k+1)
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
        if new_pctg > pctg + 9
            print(".")
            pctg = new_pctg
        end
    end # next p

    println("")
end

type WiltonInts84Strat
	outer_quad_points
end


function momintegrals!(z, op, g, f, T, τ, σ, tmax, ι, qr::WiltonInts84Strat)

    XW = qr.outer_quad_points
    for p in 1 : length(XW)
        x = XW[p].point
        w = XW[p].weight
		innerintegrals!(z, op, x, tmax, g, f, T, τ, σ, ι, qr, w)
    end # next quadrature point

end
