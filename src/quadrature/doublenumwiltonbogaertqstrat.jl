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

function integrate!(op::IntegralOperator, g::RTRefSpace, f::RTRefSpace, i, τ, j, σ, qd,
    qs::DoubleNumWiltonBogaertQStrat,
    out=nothing, test_space=nothing, tptr=nothing, trial_space=nothing, bptr=nothing;
    action::QuadRuleAction)

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

    if hits == 3
        qrule = BogaertSelfPatchStrategy(5)
        return integrate!(action, out, op, test_space, tptr, τ, trial_space, bptr, σ, qrule)
    end
    if hits == 2
        qrule = BogaertEdgePatchStrategy(8, 4)
        return integrate!(action, out, op, test_space, tptr, τ, trial_space, bptr, σ, qrule)
    end
    if hits == 1
        qrule = BogaertPointPatchStrategy(2, 3)
        return integrate!(action, out, op, test_space, tptr, τ, trial_space, bptr, σ, qrule)
    end
    rmin = xmin/k
    if xmin < xtol
        qrule = WiltonSERule(
            qd.tpoints[1,i],
            DoubleQuadRule(
                qd.tpoints[2,i],
                qd.bpoints[2,j],
            ),
        )
        return integrate!(action, out, op, test_space, tptr, τ, trial_space, bptr, σ, qrule)
    end
    qrule = DoubleQuadRule(
      qd.tpoints[1,i],
      qd.bpoints[1,j],
    )
    return integrate!(action, out, op, test_space, tptr, τ, trial_space, bptr, σ, qrule)
  end
