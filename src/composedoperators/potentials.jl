#D: dimension of manifold over which integration is performed
abstract type PotentialOperator end
struct PotentialIntegralOperator{D,U,V,W} <: PotentialOperator
    kernel::U
    op2::V
    bfunc::W
end
struct LinearCombinationPotentialOperator <: PotentialOperator
    pots
    coeffs
end
+(a::PotentialOperator,b::PotentialOperator) = LinearCombinationPotentialOperator([a,b],[1.0,1.0])
+(a::PotentialOperator,b::LinearCombinationPotentialOperator) = LinearCombinationPotentialOperator([a;b.pots],[1.0;b.coeffs])
+(a::LinearCombinationPotentialOperator,b::PotentialOperator) = LinearCombinationPotentialOperator([a.pots;b],[a.coeffs;1.0])
+(a::LinearCombinationPotentialOperator,b::LinearCombinationPotentialOperator) = LinearCombinationPotentialOperator([a.pots;b.pots],[a.coeffs;b.coeffs])
*(a::Number,b::PotentialOperator) = LinearCombinationPotentialOperator([b],[a])
-(a::PotentialOperator) = LinearCombinationPotentialOperator([a],[-1.0])
-(a::PotentialOperator,b::PotentialOperator) = a + (-1)*b

PotentialIntegralOperator{D}(a,b,c) where {D} = PotentialIntegralOperator{D,typeof(a),typeof(b),typeof(c)}(a,b,c)
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
function potential(op::PotentialIntegralOperator, points, coeffs, basis::DirectProductSpace; 
    type=SVector{3,ComplexF64},
	quadstrat=defaultquadstrat(op, basis))

    return potential(PotentialIntegralOperatorKern(op.kernel,op.op2),points,coeffs,op.bfunc(basis);type,quadstrat)
end

defaultquadstrat(op::PotentialIntegralOperator, basis) = SingleNumQStrat(6)
quaddata(op::PotentialIntegralOperatorKern,rs,els,qs::SingleNumQStrat) = quadpoints(rs,els,(qs.quad_rule,))
quadrule(op::PotentialIntegralOperatorKern,refspace,p,y,q,el,qdata,qs::SingleNumQStrat) = qdata[1,q]
kernelvals(op::PotentialIntegralOperatorKern,y,p) = nothing


function potential(op::LinearCombinationPotentialOperator, points, coeffs, basis; 
    type=SVector{3,ComplexF64})

    return sum(c*potential(pot,points,coeffs,basis;type) for (pot,c) in zip(op.pots,op.coeffs))
end

# suport for application on direct product spaces 
struct DirectProductPotentialOperatorTerm
    pot 
    hilbertvector
end
struct DirectProductPotentialOperator <: PotentialOperator
    terms
    coeffs
end
Base.getindex(a::PotentialOperator,i::Variational.HilbertVector) = DirectProductPotentialOperator([DirectProductPotentialOperatorTerm(a,i)],[1.0])
+(a::DirectProductPotentialOperator,b::DirectProductPotentialOperator) = DirectProductPotentialOperator([a.terms;b.terms],[a.coeffs;b.coeffs])
*(a::Number,b::DirectProductPotentialOperator) = DirectProductPotentialOperator(b.terms,a*b.coeffs)
-(a::DirectProductPotentialOperator) = DirectProductPotentialOperator(a.terms,-a.coeffs)

function potential(op::DirectProductPotentialOperator, points, coeffs::BlockedVector, basis::DirectProductSpace; type=SVector{3,ComplexF64})
    return sum(c*potential(op.pot, points, coeffs[op.hilbertvector],basis.factors[op.hilbertvector.idx];type) for (op,c) in zip(op.terms,op.coeffs))
end