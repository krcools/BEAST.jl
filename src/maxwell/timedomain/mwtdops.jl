

mutable struct MWSingleLayerTDIO{T} <: RetardedPotential{T}
    "speed of light in medium"
    speed_of_light::T
    "weight of the weakly singular term"
    ws_weight::T
    "weight of the hypersingular term"
    hs_weight::T
    "number of temporal differentiations in the weakly singular term"
    ws_diffs::Int
    "number of temporal differentiations in the hyper singular term"
    hs_diffs::Int
end

function Base.:*(a::Number, op::MWSingleLayerTDIO)
	@info "scalar product a * op (SL)"
	MWSingleLayerTDIO(
		op.speed_of_light,
		a * op.ws_weight,
		a * op.hs_weight,
		op.ws_diffs,
		op.hs_diffs)
end

mutable struct MWDoubleLayerTDIO{T} <: RetardedPotential{T}
    speed_of_light::T
    weight::T
    num_diffs::Int
end

function Base.:*(a::Number, op::MWDoubleLayerTDIO)
	@info "scalar product a * op (DL)"
	MWDoubleLayerTDIO(
		op.speed_of_light,
		a * op.weight,
		op.num_diffs)
end

mutable struct MWDoubleLayerTransposedTDIO{T} <: RetardedPotential{T}
	speed_of_light::T
    weight::T
    num_diffs::Int
end

function Base.:*(a::Number, op::MWDoubleLayerTransposedTDIO)
	@info "scalar product a * op (DL)"
	MWDoubleLayerTransposedTDIO(
		op.speed_of_light,
		a * op.weight,
		op.num_diffs)
end

MWSingleLayerTDIO(;speedoflight) = MWSingleLayerTDIO(speedoflight,-1/speedoflight,-speedoflight,2,0)
MWDoubleLayerTDIO(;speedoflight) = MWDoubleLayerTDIO(speedoflight, one(speedoflight), 0)


module TDMaxwell3D
import ...BEAST

function singlelayer(;speedoflight, numdiffs=0)
	@assert numdiffs >= 0
	numdiffs == 0 && return BEAST.integrate(BEAST.MWSingleLayerTDIO(speedoflight,-1/speedoflight,-speedoflight,2,0))
	return BEAST.MWSingleLayerTDIO(speedoflight,-1/speedoflight,-speedoflight,2+numdiffs-1,numdiffs-1)
end

function doublelayer(;speedoflight, numdiffs=0)
	@assert numdiffs >= -1
	numdiffs == -1 && BEAST.integrate(BEAST.MWDoubleLayerTDIO(speedoflight,1.0,0))
	return BEAST.MWDoubleLayerTDIO(speedoflight,1.0,numdiffs)
end

end # module TDMaxwell3D

export TDMaxwell3D

defaultquadstrat(::MWSingleLayerTDIO, tfs, bfs) = OuterNumInnerAnalyticQStrat(3)

function quaddata(op::MWSingleLayerTDIO, testrefs, trialrefs, timerefs,
        testels, trialels, timeels, quadstrat::OuterNumInnerAnalyticQStrat)

    dmax = numfunctions(timerefs)-1
    bn = binomial.((0:dmax),(0:dmax)')

    V = eltype(testels[1].vertices)
    ws = WiltonInts84.workspace(V)
    # quadpoints(testrefs, testels, (3,)), bn, ws
    quadpoints(testrefs, testels, (quadstrat.outer_rule,)), bn, ws
end


quadrule(op::MWSingleLayerTDIO, testrefs, trialrefs, timerefs,
        p, testel, q, trialel, r, timeel, qd, ::OuterNumInnerAnalyticQStrat) = WiltonInts84Strat(qd[1][1,p],qd[2],qd[3])


struct TransposedStorage{F}
	store::F
end

@inline (f::TransposedStorage)(v,m,n,k) = f.store(v,n,m,k)


function assemble!(dl::MWDoubleLayerTDIO, W::SpaceTimeBasis, V::SpaceTimeBasis, store,
    threading::Type{Threading{:multi}}; quadstrat=defaultquadstrat(dl,W,V))

	X, T = spatialbasis(W), temporalbasis(W)
	Y, U = spatialbasis(V), temporalbasis(V)
	if CompScienceMeshes.refines(geometry(Y), geometry(X))
		@assert !CompScienceMeshes.refines(geometry(X), geometry(Y))
		W = Y⊗T
		V = X⊗U
		op = MWDoubleLayerTransposedTDIO(dl.speed_of_light, dl.weight, dl.num_diffs)
		assemble!(op, W, V, store)
		return
	end

	P = Threads.nthreads()
	Y, S = spatialbasis(W), temporalbasis(W)
	splits = [round(Int,s) for s in range(0, stop=numfunctions(Y), length=P+1)]

	@info "Starting assembly with $P threads:"
	Threads.@threads for i in 1:P
		lo, hi = splits[i]+1, splits[i+1]
		lo <= hi || continue
		Y_p = subset(Y, lo:hi)
		store2 = (v,m,n,k) -> store(v,lo+m-1,n,k)
		assemble_chunk!(dl, Y_p ⊗ S, V, store2; quadstrat)
	end

	# return assemble_chunk!(dl, W, V, store1)
end

defaultquadstrat(::MWDoubleLayerTDIO, tfs, bfs) = OuterNumInnerAnalyticQStrat(3)

function quaddata(op::MWDoubleLayerTDIO, testrefs, trialrefs, timerefs,
        testels, trialels, timeels, quadstrat::OuterNumInnerAnalyticQStrat)

    dmax = numfunctions(timerefs)-1
    bn = binomial.((0:dmax),(0:dmax)')

    V = eltype(testels[1].vertices)
    ws = WiltonInts84.workspace(V)

    quadpoints(testrefs, testels, (quadstrat.outer_rule,)), bn, ws
end

quadrule(op::MWDoubleLayerTDIO, testrefs, trialrefs, timerefs,
    p, testel, q, trialel, r, timeel, qd, quadstrat::OuterNumInnerAnalyticQStrat) =
        WiltonInts84Strat(qd[1][1,p],qd[2],qd[3])


defaultquadstrat(::MWDoubleLayerTransposedTDIO, tfs, bfs) = OuterNumInnerAnalyticQStrat(3)

function quaddata(op::MWDoubleLayerTransposedTDIO,
		testrefs, trialrefs, timerefs,
        testels, trialels, timeels, quadstrat::OuterNumInnerAnalyticQStrat)

    dmax = numfunctions(timerefs)-1
    bn = binomial.((0:dmax),(0:dmax)')

    V = eltype(testels[1].vertices)
    ws = WiltonInts84.workspace(V)
    quadpoints(testrefs, testels, (quadstrat.outer_rule,)), bn, ws
end

quadrule(op::MWDoubleLayerTransposedTDIO, testrefs, trialrefs, timerefs,
    p, testel, q, trialel, r, timeel, qd, quadstrat::OuterNumInnerAnalyticQStrat) =
        WiltonInts84Strat(qd[1][1,p],qd[2],qd[3])

function momintegrals!(z, op::MWDoubleLayerTransposedTDIO,
	g, f, T, τ, σ, ι, qr::WiltonInts84Strat)

	op1 = MWDoubleLayerTDIO(op.speed_of_light, op.weight, op.num_diffs)
	momintegrals!(z, op1, g, f, T, τ, σ, ι, qr::WiltonInts84Strat)
	w = similar(z)
	permutedims!(w, z, [2,1,3])

end

@inline function tmRoR(d, iG, bns)
    sgn = isodd(d) ? -1 : 1
    r = sgn * iG[d+2]
end

# build
# ``\int (D-R)^d/R (y-b) dy`` from
# ``(ξ-b) \int R^k dy`` and
# ``\int R^k (y-ξ) dy``
@inline function tmRoRf(d, iG, iGξy, bξ, bns)
    sgn = isodd(d) ? -1 : 1
    iGf = iGξy[d+2] + bξ * iG[d+2]
    r = sgn * iGf
end

"""
    Q = qd(T,dh,::Val{N})

Q[k] is the factor in front resulting from differentiating t^(k-1) dh times.
"""
@generated function qh(::Type{T},dh,n::Type{Val{N}}) where {N,T}
    xp = quote end
    for k in 1:N
        qk = Symbol(:q,k)
        d = k-1
        xp1 = quote
            $(qk) = one($T)
            for v in 0 : dh-1
                $(qk) *= ($d-v)
            end
        end
        append!(xp.args, xp1.args)
    end
    xp1 = :(())
    for k in 1:N
        qk = Symbol(:q,k)
        push!(xp1.args, :($qk))
    end
    push!(xp.args, xp1)
    return xp
end


function innerintegrals!(zl, op::MWSingleLayerTDIO,
        p, # test_point, test_time
        U, V, W, # local_test_space, local_trial_space, local_temporal_space
        τ, σ, ι, # test_element, trial_element, spherial_shell
        qr, w)   # inner_quadrature_rule, outer_quadrature_weight

	T = typeof(w)

    sol = op.speed_of_light
    #Rmax = sol * tmax

    dx = w
    x = cartesian(p)
    n = cross(σ[1]-σ[3],σ[2]-σ[3])
    n /= norm(n)
    ξ = x - ((x-σ[1]) ⋅ n) * n

    r = ι[1]
    R = ι[2]

    @assert r < R
    @assert degree(W) <= 3

    ∫G, ∫Gξy, = WiltonInts84.wiltonints(σ[1],σ[2],σ[3],x,r,R,Val{2},qr.workspace)

    αg = 1 / volume(τ) / 2
	αf = 1 / volume(σ) / 2
	αG = 1 / 4π
	α = αg * αf * αG * op.ws_weight * dx
	β = 4 * αg * αf * αG * op.hs_weight * dx

	ds = op.ws_diffs
	dh = op.hs_diffs

    qhs = qh(T,dh,Val{4})
    qss = qh(T,ds,Val{4})

    bn = qr.binomials

    #solpowers = collect(sol^p for p ∈ 0:numfunctions(W)-1)
    sol2 = sol*sol
    sol3 = sol2*sol
    sol4 = sol3*sol
    sol5 = sol4*sol
    solpowers = (one(sol), sol, sol2, sol3, sol4, sol5)

    udim = numfunctions(U, domain(τ))
    vdim = numfunctions(V, domain(σ))

    for i in 1 : udim
        a = τ[i]
        g = (x-a)
        for j in 1 : vdim
            b = σ[j]; bξ = ξ-b
            for k in 1 : numfunctions(W)
				d = k-1 # ranges from 0 to numfunctions(W)-1
				sgn = isodd(d) ? -1 : 1
				# hyper singular contribution
				if d >= dh
                    @assert dh == 0
                    q = qhs[k]
                    Ih = tmRoR(d-dh, ∫G, bn) # \int (cTmax-R)^(d-dh)/R dy
                    #zl[i,j,k] += β * q * Ih / sol^(d-dh)
                    zl[i,j,k] += β * q * Ih / solpowers[d-dh+1]
				end
				# weakly singular contribution
				if d >= ds
                    q = qss[k]
                    Is = tmRoRf(d-ds, ∫G, ∫Gξy, bξ, bn) # \int (cTmax-R)^(d-ds)/R (y-b) dy
                    #zl[i,j,k] += α * q * (g ⋅ Is) / sol^(d-ds)
                    zl[i,j,k] += α * q * (g ⋅ Is) / solpowers[d-ds+1]
				end
            end
        end
    end

end



function innerintegrals!(z, op::MWDoubleLayerTDIO,
    p,
    U, V, W,
    τ, σ, ι,
    qr, w)

	T = typeof(w)

    sol = op.speed_of_light

    dx = w
    x = cartesian(p)
    n = cross(σ[1]-σ[3],σ[2]-σ[3])
    n /= norm(n)
    ξ = x - ((x-σ[1]) ⋅ n) * n

    r = ι[1]
    R = ι[2]

    @assert r < R
    @assert degree(W) <= 3

    #N = max(degree(W), 0)
    ∫G, ∫Gξy, ∫∇G = WiltonInts84.wiltonints(σ[1],σ[2],σ[3],x,r,R,Val{2},qr.workspace)
	@assert isapprox(∫∇G[2] , point(0,0,0), atol=1e-8)

    αg = 1 / volume(τ) / 2
	αf = 1 / volume(σ) / 2
	αG = 1 / 4 / π
	α = αG * op.weight * dx # * αg * αf

	ds = op.num_diffs

    @inline function tmGRoR(d, iGG)
        sgn = isodd(d) ? -1 : 1
        r = sgn * iGG[d+1]
    end

    Ux = U(p)
    Vx = αf * @SVector[(x-σ[1]), (x-σ[2]), (x-σ[3])]

    udim = numfunctions(U, domain(τ))
    vdim = numfunctions(V, domain(σ))

    for i in 1 : udim
        # a = τ[i]
        # g = αg * (x-τ[i])
        g = Ux[i].value
        for j in 1 : vdim
            # b = σ[j]
            # f = αf * (x-σ[j])
            # f = Vx[j].value
            f = Vx[j]
            fxg = f × g
            for k in 1 : numfunctions(W)
				d = k-1
				sgn = isodd(d) ? -1 : 1
				if d >= ds
					q = one(T)
					for p in 0 : ds-1
						q *= (d-p)
					end
                    # @assert q == 1
                    z[i,j,k] += -α * q * ( fxg ⋅ tmGRoR(d-ds, ∫∇G) ) / sol^(d-ds)
				end
            end
        end
    end

end
