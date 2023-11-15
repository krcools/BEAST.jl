abstract type FarField end

defaultquadstrat(op::FarField, basis) = SingleNumQStrat(3)
defaultquadstrat(op::FarField, basis::DirectProductSpace) = defaultquadstrat(op, basis.factors[1])

quaddata(op::FarField,rs,els,qs::SingleNumQStrat) = quadpoints(rs,els,(qs.quad_rule,))
quadrule(op::FarField,refspace,p,y,q,el,qdata,qs::SingleNumQStrat) = qdata[1,q]