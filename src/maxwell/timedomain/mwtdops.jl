

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

mutable struct MWDoubleLayerTDIO{T} <: RetardedPotential{T}
    speed_of_light::T
    weight::T
    num_diffs::Int
end

MWSingleLayerTDIO(;speedoflight) = MWSingleLayerTDIO(speedoflight,-1/speedoflight,-speedoflight,2,0)
MWDoubleLayerTDIO(;speedoflight) = MWDoubleLayerTDIO(speedoflight, one(speedoflight), 0)

module TimeDomain
module Maxwell3D
import ...BEAST
SingleLayer(;speedoflight) = BEAST.MWSingleLayerTDIO(speedoflight,-1/speedoflight,-speedoflight,2,0)
end
end

function quaddata(op::MWSingleLayerTDIO, testrefs, trialrefs, timerefs,
        testels, trialels, timeels)

    dmax = numfunctions(timerefs)-1
    bn = binomial.((0:dmax),(0:dmax)')

    V = eltype(testels[1].vertices)
    ws = WiltonInts84.workspace(V)
    quadpoints(testrefs, testels, (3,)), bn, ws

end


quadrule(op::MWSingleLayerTDIO, testrefs, trialrefs, timerefs,
        p, testel, q, trialel, r, timeel, qd) = WiltonInts84Strat(qd[1][1,p],qd[2],qd[3])

function quaddata(op::MWDoubleLayerTDIO, testrefs, trialrefs, timerefs,
        testels, trialels, timeels)

    dmax = numfunctions(timerefs)-1
    bn = binomial.((0:dmax),(0:dmax)')

    V = eltype(testels[1].vertices)
    ws = WiltonInts84.workspace(V)
    quadpoints(testrefs, testels, (3,)), bn, ws
end

quadrule(op::MWDoubleLayerTDIO, testrefs, trialrefs, timerefs,
        p, testel, q, trialel, r, timeel, qd) = WiltonInts84Strat(qd[1][1,p],qd[2],qd[3])

@inline function tmRoR(d, iG, bns)
    # r = zero(t)
    # for q in 0:d
    #     sgn = isodd(q) ? -1 : 1
    #     bn = bns[d+1,q+1]
    #     #r += binomial(d,q) * sgn * t^(d-q) * iG[q+2]
    #     r += bn * sgn * t^(d-q) * iG[q+2]
    # end
    # r
    sgn = isodd(d) ? -1 : 1
    r = sgn * iG[d+2]
end

# build
# ``\int (D-R)^d/R (y-b) dy`` from
# ``(ξ-b) \int R^k dy`` and
# ``\int R^k (y-ξ) dy``
@inline function tmRoRf(d, iG, iGξy, bξ, bns)
    # r = zero(bξ)
    # for q in 0:d
    #     sgn = isodd(q) ? -1 : 1
    #     i = q+2
    #     iGf = iGξy[i] + bξ * iG[i]
    #     bn = bns[d+1,q+1]
    #     #r += binomial(d,q) * sgn * t^(d-q) * iGf
    #     r += bn * sgn * t^(d-q) * iGf
    # end
    # r
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

    for i in 1 : numfunctions(U)
        a = τ[i]
        g = (x-a)
        for j in 1 : numfunctions(V)
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
    #Rmax = sol * tmax

    dx = w
    x = cartesian(p)
    #n = normal(σ)
    n = cross(σ[1]-σ[3],σ[2]-σ[3])
    n /= norm(n)
    ξ = x - ((x-σ[1]) ⋅ n) * n

    r = ι[1]
    R = ι[2]

    @assert r < R
    @assert degree(W) <= 3

    #N = max(degree(W), 0)
    ∫G, ∫Gξy, ∫∇G = WiltonInts84.wiltonints(σ[1],σ[2],σ[3],x,r,R,Val{2},qr.workspace)

    αg = 1 / volume(τ) / 2
	αf = 1 / volume(σ) / 2
	αG = 1 / 4 / π
	α = αg * αf * αG * op.weight * dx

	ds = op.num_diffs

    @inline function tmRoR(d, iGG)
        # r = zero(t)
        # for q in 0:d
        #     sgn = isodd(q) ? -1 : 1
        #     r += binomial(d,q) * sgn * t^(d-q) * iGG[q+1]
        # end
        # r
        sgn = isodd(d) ? -1 : 1
        r = sgn * iGG[d+1]
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
                    z[i,j,k] += -α * q * ( fxg ⋅ tmRoR(d-ds, ∫∇G) ) / sol^(d-ds)
				end
            end
        end
    end

end
