

"""
Operator to compute the far field of a current distribution in the time domain.
In particular, given the current distribution ``j`` this operator allows for the computation of

```math
R =  ̂x ⋅ y
ffd = n × ∫_Γ j(r', t - R/c} dy
```

where ``̂x`` is the unit vector in the direction of observation.
Note that the assembly routing expects the observation directions to be normalised by the caller.
"""
mutable struct MWFarField3DTD{K}
  sol::K
end

function quaddata(op::MWFarField3DTD, trialrefs, timerefs, trialels, timeels)

    trialqd = quadpoints(trialrefs, trialels, (3,))

    trialqd, nothing

end

function quadrule(op::MWFarField3DTD,trialrefs, timerefs,
        p, testel, q, trialel, r, timeel, qd)
    qd[1][1,q], nothing
end

kernelvals(op::MWFarField3DTD,test,source) = dot(test,cartesian(source))/op.sol

function integrand(op::MWFarField3DTD,krn,testels, trialvals, t, T)
    τ = krn
    timevals = T(t - τ)

    testels × (trialvals[1] * timevals)
end
