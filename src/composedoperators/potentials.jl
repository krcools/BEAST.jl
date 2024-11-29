#D: dimension of manifold over which integration is performed
struct PotentialIntegralOperator{D}
    kernel
    op2
    bfunc
end

struct PotentialIntegralOperatorKern{U,V}
    kernel::U
    op2::V
end

function integrand(a::PotentialIntegralOperatorKern,y,krn,f,x)
    return a.op2(a.kernel(y,x),f)
end
function potential(op::PotentialIntegralOperator, points, coeffs, basis; 
    type=SVector{3,ComplexF64},
	quadstrat=defaultquadstrat(op, basis))

    return potential(PotentialIntegralOperatorKern(op.kernel,op.op2),points,coeffs,op.bfunc(basis);type,quadstrat)
end

