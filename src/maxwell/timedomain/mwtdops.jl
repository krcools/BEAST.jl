export MWSingleLayerTDIO
export MWDoubleLayerTDIO

type MWSingleLayerTDIO{T} <: RetardedPotential{T}
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

type MWDoubleLayerTDIO{T} <: RetardedPotential{T}
    speed_of_light::T
    weight::T
    num_diffs::Int
end

MWSingleLayerTDIO(c) = MWSingleLayerTDIO(c,-1/c,-c,2,0)
MWDoubleLayerTDIO(c) = MWDoubleLayerTDIO(c, one(c), 0)

quaddata(op::MWSingleLayerTDIO, testrefs, trialrefs, timerefs,
    testels, trialels, timeels) = quadpoints(testrefs, testels, (3,))
quadrule(op::MWSingleLayerTDIO, testrefs, trialrefs, timerefs,
        p, testel, q, trialel, r, timeel, qd) = WiltonInts84Strat(qd[1,p])

quaddata(op::MWDoubleLayerTDIO, testrefs, trialrefs, timerefs,
        testels, trialels, timeels) = quadpoints(testrefs, testels, (3,))
quadrule(op::MWDoubleLayerTDIO, testrefs, trialrefs, timerefs,
        p, testel, q, trialel, r, timeel, qd) = WiltonInts84Strat(qd[1,p])


function innerintegrals!(zl, op::MWSingleLayerTDIO,
    p, tmax, # test_point, test_time
    U, V, W, # local_test_space, local_trial_space, local_temporal_space
    τ, σ, ι, # test_element, trial_element, spherial_shell
    qr, w)   # inner_quadrature_rule, outer_quadrature_weight

	T = typeof(w)

    dx = w
    x = cartesian(p)
    #n = normal(σ)
    n = cross(σ[1]-σ[3],σ[2]-σ[3])
    n /= norm(n)
    ξ = x - ((x-σ[1]) ⋅ n) * n

    r = ι[1]
    R = ι[2]

    @assert r < R

    N = max(degree(W), 0)
    ∫G, ∫Gξy, = WiltonInts84.wiltonints(σ[1],σ[2],σ[3],x,r,R,Val{N-1})

    αg = 1 / volume(τ) / 2
	αf = 1 / volume(σ) / 2
	αG = 1 / 4π
	α = αg * αf * αG * op.ws_weight * dx
	β = 4 * αg * αf * αG * op.hs_weight * dx

	ds = op.ws_diffs
	dh = op.hs_diffs

    @inline function tmRoR(d, t, iG)
        r = zero(t)
        for q in 0:d
            sgn = isodd(q) ? -1 : 1
            r += binomial(d,q) * sgn * t^(d-q) * iG[q+2]
        end
        r
    end

    @inline function tmRoRf(d, t, iG, iGξy, bξ)
        r = zero(t)
        for q in 0:d
            sgn = isodd(q) ? -1 : 1
            i = q+2
            iGf = iGξy[i] + bξ * iG[i]
            r += binomial(d,q) * sgn * t^(d-q) * iGf
        end
        r
    end

    for i in 1 : numfunctions(U)
        a = τ[i]
        g = (x-a)
        for j in 1 : numfunctions(V)
            b = σ[j]; bξ = ξ-b
            #∫Gf = ∫Gξy + (ξ-b) * ∫G
            #∫Gf = [∫Gξy[i] + bξ*∫G[i] for i in eachindex(∫G)]
            for k in 1 : numfunctions(W)
				d = k-1
				sgn = isodd(d) ? -1 : 1
				# hyper singular contribution
				if d >= dh
					q = one(T)
					for p in 0 : dh-1
						q *= (d-p)
					end
                    @assert dh == 0
                    zl[i,j,k] += β * q * tmRoR(d-dh, tmax, ∫G)
					#zl[i,j,k] += β * sgn * q * ∫G[d+2-dh]
				end
				# weakly singular contribution
				if d >= ds
					q = one(T)
					for p in 0 : ds-1
						q *= (d-p)
					end
					#∫Gf = ∫Gξy[d+2-ds] + (ξ-b)*∫G[d+2-ds]
                    zl[i,j,k] += α * q * (g ⋅ tmRoRf(d-ds, tmax, ∫G, ∫Gξy, bξ))
					#zl[i,j,k] += α * sgn * q * (g ⋅ ∫Gf)
				end
            end
        end
    end

end



function innerintegrals!(z, op::MWDoubleLayerTDIO,
    p, tmax,
    U, V, W,
    τ, σ, ι,
    qr, w)

	T = typeof(w)

    dx = w
    x = cartesian(p)
    #n = normal(σ)
    n = cross(σ[1]-σ[3],σ[2]-σ[3])
    n /= norm(n)
    ξ = x - ((x-σ[1]) ⋅ n) * n

    r = ι[1]
    R = ι[2]

    @assert r < R

    N = max(degree(W), 0)
    ∫G, ∫Gξy, ∫∇G = WiltonInts84.wiltonints(σ[1],σ[2],σ[3],x,r,R,Val{N-1})

    αg = 1 / volume(τ) / 2
	αf = 1 / volume(σ) / 2
	αG = 1 / 4 / π
	α = αg * αf * αG * op.weight * dx

	ds = op.num_diffs

    @inline function tmRoR(d, t, iGG)
        r = zero(t)
        for q in 0:d
            sgn = isodd(q) ? -1 : 1
            #r += binomial(d,q) * sgn * t^(d-q) * iG[q+2]
            r += binomial(d,q) * sgn * t^(d-q) * iGG[q+1]
        end
        r
    end

    for i in 1 : numfunctions(U)
        a = τ[i]
        g = (x-a)
        for j in 1 : numfunctions(V)
            b = σ[j]
            f = (x-b)
            fxg = f × g
            for k in 1 : numfunctions(W)
				d = k-1
				sgn = isodd(d) ? -1 : 1
				if d >= ds
					q = one(T)
					for p in 0 : ds-1
						q *= (d-p)
					end
                    @assert q == 1
                    #z[i,j,k] += -α * sgn * q * (f × g) ⋅ ∫∇G[d+1-ds]
                    z[i,j,k] += -α * q * ( fxg ⋅ tmRoR(d-ds, tmax, ∫∇G) )
                    #d == 1 && @assert z[i,j,k] == 0
				end
            end
        end
    end

end
