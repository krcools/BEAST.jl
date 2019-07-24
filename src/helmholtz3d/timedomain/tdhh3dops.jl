abstract type HH3DTDBIO{T} <: RetardedPotential{T} end

struct HH3DSingleLayerTDBIO{T} <: HH3DTDBIO{T}
    "speed of light"
    speed_of_light::T
    "number of temporal differentiations"
    num_diffs::Int
end
HH3DSingleLayerTDBIO(c) = HH3DSingleLayerTDBIO(c, 0)

struct HH3DHyperSingularTDBIO{T} <: HH3DTDBIO{T}
    speed_of_light::T
    weight_of_weakly_singular_term::T
    weight_of_hyper_singular_term::T
    num_diffs_on_weakly_singular_term::Int
    num_diffs_on_hyper_singular_term::Int
end
function HH3DHyperSingularTDBIO(;speed_of_light, numdiffs)
    id = one(speed_of_light)
    HH3DHyperSingularTDBIO(speed_of_light, id, id, numdiffs+1, numdiffs-1)
end

# See: ?BEAST.quaddata for help
function quaddata(operator::HH3DTDBIO,
        test_local_space, trial_local_space, time_local_space,
        test_element, trial_element, time_element)

    dmax = numfunctions(time_local_space)-1
    bn = binomial.((0:dmax),(0:dmax)')

    V = eltype(test_element[1].vertices)
    ws = WiltonInts84.workspace(V)
    quadpoints(test_local_space, test_element, (3,)), bn, ws

end


# See: ?BEAST.quadrule for help
function quadrule(operator::HH3DTDBIO,
        test_local_space, trial_local_space, time_local_space,
        p, test_element, q, trial_element, r, time_element,
        quad_data)

    # WiltonInts84Strat(quad_data[1,p])
    qd = quad_data
    WiltonInts84Strat(qd[1][1,p],qd[2],qd[3])

end


function innerintegrals!(zlocal, operator::HH3DSingleLayerTDBIO,
        test_point,
        test_local_space, trial_local_space, time_local_space,
        test_element, trial_element, time_element,
        quad_rule, quad_weight)

    dx = quad_weight
    x = cartesian(test_point)
    n = normal(test_point)

    a = trial_element[1]
    ξ = x - dot(x -a, n) * n

    r = time_element[1]
    R = time_element[2]
    @assert r < R

    N = max(degree(time_local_space), 1)
    ∫G, = WiltonInts84.wiltonints(
        trial_element[1],
        trial_element[2],
        trial_element[3],
        x, r, R, Val{N-1},quad_rule.workspace)

    a = dx / (4*pi)
    D = operator.num_diffs
    @assert D == 0
    @assert numfunctions(test_local_space)  == 1
    @assert numfunctions(trial_local_space) == 1

    # # aux fn to compute \int (t-R)^d / R
    # @inline function tmRoR(d, t, iG)
    #     r = zero(t)
    #     for q in 0:d
    #         sgn = isodd(q) ? -1 : 1
    #         r += binomial(d,q) * sgn * t^(d-q) * iG[q+2]
    #     end
    #     r
    # end

    @inline function tmRoR(d, iG, bns)
        sgn = isodd(d) ? -1 : 1
        r = sgn * iG[d+2]
    end

    bns = quad_rule.binomials

    for k in 1 : numfunctions(time_local_space)
        d = k - 1
        if d >= D
            q = 1
            for p in 0 : D-1
                q *= (d-p)
            end
            # zlocal[1,1,k] += a * q * tmRoR(d-D, test_time, ∫G)
            zlocal[1,1,k] += a * q * tmRoR(d-D, ∫G, bns)
            # sgn = isodd(d) ? -1 : 1
            # zlocal[1,1,k] += a * sgn * q * ∫G[d+2-D]
        end
    end # k
    # comment

end


function innerintegrals!(zlocal, operator::HH3DHyperSingularTDBIO,
        test_point,
        test_local_space, trial_local_space, time_local_space,
        test_element, trial_element, time_element,
        quad_rule, quad_weight)

    dx = quad_weight
    x = cartesian(test_point)
    # n = normal(test_point)

    nx = normal(test_point)
    ny = normal(trial_element)
    ndotn = dot(nx,ny)

    a = trial_element[1]
    ξ = x - dot(x -a, ny) * ny

    r = time_element[1]
    R = time_element[2]
    @assert r < R

    N = max(degree(time_local_space), 1)
    ∫G, ∫Gξy, = WiltonInts84.wiltonints(
        trial_element[1],
        trial_element[2],
        trial_element[3],
        x, r, R, Val{N-1},quad_rule.workspace)

    # a = dx / (4*pi)
    # D = operator.num_diffs
    # @assert D == 0
    @assert numfunctions(test_local_space)  == 3
    @assert numfunctions(trial_local_space) == 3

    @inline function tmRoR(d, iG)
        r = (isodd(d) ? -1 : 1) * iG[d+2]
    end

    @inline function tmRoRf(d, ∫G, ∫Gξy, bξ, h, m)
        (isodd(d) ? -1 : 1) * ((1 - h*dot(m,bξ)) * ∫G[d+2] - h*dot(m, ∫Gξy[d+2]))
    end

    test_values = test_local_space(test_point, Val{:withcurl})
    trial_values = trial_local_space(center(trial_element), Val{:withcurl})

    # weakly singular term
    α = dx / (4π) * operator.weight_of_weakly_singular_term
    Ds = operator.num_diffs_on_weakly_singular_term
    for i in 1 : numfunctions(test_local_space)
        g, curlg = test_values[i]
        for j in 1 : numfunctions(trial_local_space)
            b = trial_element[j]
            opp_edge = trial_element[mod1(j+2,3)] - trial_element[mod1(j+1,3)]
            h = norm(opp_edge)/2/volume(trial_element)
            m = normalize(cross(opp_edge, ny))
            for k in 1 : numfunctions(time_local_space)
                d = k-1
                d < Ds && continue
                q = reduce(*, d-Ds+1:d ,init=1)
                zlocal[i,j,k] += α * q * g * ndotn * tmRoRf(d-Ds, ∫G, ∫Gξy, ξ-b, h, m)
            end
        end
    end

    # Hyper-singular term
    β = dx / (4π) * operator.weight_of_hyper_singular_term
    Dh = operator.num_diffs_on_hyper_singular_term
    for i in 1 : numfunctions(test_local_space)
        g, curlg = test_values[i]
        for j in 1 : numfunctions(trial_local_space)
            _, curlf = trial_values[j]
            for k in 1 : numfunctions(time_local_space)
                d = k - 1
                d < Dh && continue
                q = reduce(*, d-Dh+1:d ,init=1)
                zlocal[i,j,k] += β * q * dot(curlg, curlf) * tmRoR(d-Dh, ∫G)
            end
        end
    end
end
