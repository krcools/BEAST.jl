struct DoubleNumWiltonBogaertQStrat{R} <: AbstractQuadStrat
    outer_rule_far::R
    inner_rule_far::R
    outer_rule_near::R
    inner_rule_near::R
end

function quaddata(op::IntegralOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts, trial_charts, qs::DoubleNumWiltonBogaertQStrat)

    T = coordtype(test_charts[1])

    tqd = quadpoints(test_local_space,  test_charts,  (qs.outer_rule_far,qs.outer_rule_near))
    bqd = quadpoints(trial_local_space, trial_charts, (qs.inner_rule_far,qs.inner_rule_near))

    return (tpoints=tqd, bpoints=bqd)
end

function quadrule(op::IntegralOperator, g::RTRefSpace, f::RTRefSpace, i, τ, j, σ, qd,
    qs::DoubleNumWiltonBogaertQStrat)

    dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))
    xtol = 0.2
  
    k = norm(gamma(op))
  
    hits = 0
    xmin = xtol
    for t in τ.vertices
      for s in σ.vertices
        d = norm(t-s)
        xmin = min(xmin, k*d)
        if d < dtol
          hits +=1
          break
        end
      end
    end
  
    hits == 3   && return BogaertSelfPatchStrategy(5)
    hits == 2   && return BogaertEdgePatchStrategy(8, 4)
    hits == 1   && return BogaertPointPatchStrategy(2, 3)
    rmin = xmin/k
    xmin < xtol && return WiltonSERule(
      qd.tpoints[1,i],
      DoubleQuadRule(
        qd.tpoints[2,i],
        qd.bpoints[2,j],
      ),
    )
    return DoubleQuadRule(
      qd.tpoints[1,i],
      qd.bpoints[1,j],
    )
  
  end