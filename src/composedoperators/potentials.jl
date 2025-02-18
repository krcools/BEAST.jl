#D: dimension of manifold over which integration is performed
abstract type PotentialOperator end
struct PotentialIntegralOperator{D,U,V,W} <: PotentialOperator
    kernel::U
    op2::V
    bfunc::W
end

struct PotentialIntegralOperatorKern{U,V}
    kernel::U
    op2::V
end

function integrand(a::PotentialIntegralOperatorKern,krn,y,f,x)
    return a.op2(a.kernel(y,x),f.value)
end
function potential(op::PotentialIntegralOperator, points, coeffs, basis; 
    type=SVector{3,ComplexF64},
	quadstrat=defaultquadstrat(op, basis))

    return potential(PotentialIntegralOperatorKern(op.kernel,op.op2),points,coeffs,op.bfunc(basis);type,quadstrat)
end

defaultquadstrat(op::PotentialIntegralOperator, basis) = SingleNumQStrat(6)
quaddata(op::PotentialIntegralOperatorKern,rs,els,qs::SingleNumQStrat) = quadpoints(rs,els,(qs.quad_rule,))
quadrule(op::PotentialIntegralOperatorKern,refspace,p,y,q,el,qdata,qs::SingleNumQStrat) = qdata[1,q]
kernelvals(op::PotentialIntegralOperatorKern,y,p) = nothing