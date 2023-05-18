# This file contains kernel independent implementations of quaddata for
# the various quadrature strategies defined in quadstrats.jl

function quaddata(op::IntegralOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts, trial_charts, qs::DoubleNumSauterQstrat)

    T = coordtype(test_charts[1])

    tqd = quadpoints(test_local_space,  test_charts,  (qs.outer_rule,))
    bqd = quadpoints(trial_local_space, trial_charts, (qs.inner_rule,))
     
    leg = (
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_vert,0,1)),
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_edge,0,1)),
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_face,0,1)),
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_tetr,0,1)),
      )

    return (tpoints=tqd, bpoints=bqd, gausslegendre=leg)
end