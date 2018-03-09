
# T: fieldtype
struct HH3DSingleLayerTDBIO{T} <: RetardedPotential{T}
    "speed of light"
    speed_of_light::T
    "number of temporal differentiations"
    num_diffs
end

HH3DSingleLayerTDBIO(c) = HH3DSingleLayerTDBIO(c, 0)

# See: ?BEAST.quaddata for help
function quaddata(operator::HH3DSingleLayerTDBIO,
        test_local_space, trial_local_space, time_local_space,
        test_element, trial_element, time_element)

    quaddata(test_local_space, test_element, (3,))

end


# See: ?BEAST.quadrule for help
function quadrule(operator::HH3DSingleLayerTDBIO,
        test_local_space, trial_local_space, time_local_space,
        p, test_element, q, trial_element, r, time_element,
        quad_data)

    WiltonInts84Strat(quad_data[1,p])

end


function innerintegrals!(zlocal, operator::HH3DSingleLayerTDBIO,
        test_point, test_time,
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
        x, r, R, Val{N-1})

    a = dx / (4*pi)
    D = operator.num_diffs
    @assert D == 0
    @assert numfunctions(test_local_space)  == 1
    @assert numfunctions(trial_local_space) == 1

    # aux fn to compute \int (t-R)^d / R
    @inline function tmRoR(d, t, iG)
        r = zero(t)
        for q in 0:d
            sgn = isodd(q) ? -1 : 1
            r += binomial(d,q) * sgn * t^(d-q) * iG[q+2]
        end
        r
    end

    for k in 1 : numfunctions(time_local_space)
        d = k - 1
        if d >= D
            q = 1
            for p in 0 : D-1
                q *= (d-p)
            end
            zlocal[1,1,k] += a * q * tmRoR(d-D, test_time, ∫G)
            # sgn = isodd(d) ? -1 : 1
            # zlocal[1,1,k] += a * sgn * q * ∫G[d+2-D]
        end
    end # k

end
