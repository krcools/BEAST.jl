struct DisplacementVector{T} <: Kernel{T} 
    interior::Bool
end


function (u::DisplacementVector)(x,disp)
    if disp == 0
        return (u.interior - !u.interior) * normal(x)
    else
        return sign(disp)*normal(x)
    end
end


"""
    _trace(kernel,interior,Val{N})

    computes the jump contribution of a kernel, where N is the dimension of the potential integral.
"""
_trace(::Kernel,interior,_) = 0, 0
function _trace(::HH3DGradGreen{T},interior::Bool,::Val{2}) where {T}
    return DisplacementVector{real(T)}(interior) , real(T)(0.5)
end

# """
#     ttrace(Potential,interior::Bool;testfunction_tangential=false)

# Compute the tangential trace, -n×n× Potential, of a potential operator mapping it to a boundary operator.

# This function assumes the normalvector on the mesh to point outwards. 
# Global multi-trace structures are supported when this function is called.
# """
# function ttrace(a::PotentialIntegralOperator{D},interior::Bool,testfunctionmap = x->(n× x)×n) where {D}
#     t,f = _trace(a.kernel,interior,Val{D}())
#     doubleint = CompDoubleInt(B->testfunctionmap(B),(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
#     t!=0 && (singleint = f*CompSingleInt(B->testfunctionmap(B),(x,y)->transpose(x)*y,t,a.op2,a.bfunc))
#     if t==0
#         return doubleint
#     else
#         #eventualy a scan to be non zero can be added here.
#         return doubleint + singleint
#     end
# end
# """
#     ntrace(Potential,interior::Bool)

# Compute the normal trace, n⋅ Potential, of a potential operator mapping it to a boundary operator.

# This function assumes the normalvector on the mesh to point outwards. 
# Global multi-trace structures are supported when this function is called.
# """
# function ntrace(a::PotentialIntegralOperator{D},interior::Bool,testfunctionmap = x->x*n) where {D}
#     t,f = _trace(a.kernel,interior,Val{D}())
#     doubleint = CompDoubleInt(B-> testfunctionmap(B),(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
#     t!=0 && (singleint = f*CompSingleInt(B-> testfunctionmap(B),(x,y)->transpose(x)*y,t,a.op2,a.bfunc))
    
#     if t==0
#         return doubleint
#     else
#         #eventualy a scan to be non zero can be added here.
#         return doubleint + singleint
#     end
# end
# """
#     trace(Potential,interior::Bool)

# Compute the trace of a potential operator mapping it to a boundary operator.

# This function assumes the normalvector on the mesh to point outwards. 
# Global multi-trace structures are supported when this function is called.
# """
# function trace(a::PotentialIntegralOperator{D},interior::Bool,testfunctionmap = x->x) where {D}
#     t,f = _trace(a.kernel,interior,Val{D}())
#     doubleint = CompDoubleInt(B->testfunctionmap(B),*,a.kernel,a.op2,a.bfunc)
#     t!=0 && (singleint = f*CompSingleInt(B->testfunctionmap(B),*,t,a.op2,a.bfunc))
    
#     if t==0
#         return doubleint
#     else
#         #eventualy a scan to be non zero can be added here.
#         return doubleint + singleint
#     end
# end
# """
#     rtrace(Potential,interior::Bool)

# Compute the rotated trace, n × Potential, of a potential operator mapping it to a boundary operator.

# This function assumes the normalvector on the mesh to point outwards. 
# Global multi-trace structures are supported when this function is called.
# """
# function rtrace(a::PotentialIntegralOperator{D},interior::Bool,testfunctionmap = x->x×n) where {D}
#     t,f = _trace(a.kernel,interior,Val{D}())
#     doubleint = CompDoubleInt(B-> testfunctionmap(B),(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
#     t!=0 && (singleint = f*CompSingleInt(B-> testfunctionmap(B),(x,y)->transpose(x)*y,t,a.op2,a.bfunc))
    
#     if t==0
#         return doubleint
#     else
#         #eventualy a scan to be non zero can be added here.
#         return doubleint + singleint
#     end
# end
# """
#     pvttrace(Potential,interior::Bool;testfunction_tangential=false)

# Compute the principal value of the tangential trace, -n×n× Potential.

# A small acceleration can be gained by setting testfunction_tangential to true.
# """
# # the principal value of the trace.
# function pvttrace(a::PotentialIntegralOperator{D},testfunctionmap = x->n×(x×n);testfunction_tangential=false) where {D}
#     if testfunction_tangential
#         doubleint = CompDoubleInt(B->testfunctionmap(B),(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
#     else
#         doubleint = CompDoubleInt(B-> testfunctionmap(B),(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
#     end
#     return doubleint
# end

# """
#     pvntrace(Potential,interior::Bool)

# Compute the principal value of the normal trace, n⋅ Potential.

# This function assumes the normalvector on the mesh to point outwards. 
# """
# function pvntrace(a::PotentialIntegralOperator{D}, testfunctionmap = x->n*x) where {D}

#     doubleint = CompDoubleInt(B-> testfunctionmap(B),(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)

#     return doubleint
# end

# """
#     pvtrace(Potential,interior::Bool)

# Compute the principal value of the  trace.
# """
# function pvtrace(a::PotentialIntegralOperator{D},testfunctionmap = x->x) where {D}
#     doubleint = CompDoubleInt(B-> testfunctionmap(B),*,a.kernel,a.op2,a.bfunc)

#     return doubleint
# end

# """
#     pvrtrace(Potential,interior::Bool)

# Compute the principal value of the rotated trace, n× Potential.

# This function assumes the normalvector on the mesh to point outwards. 
# """
# function pvrtrace(a::PotentialIntegralOperator{D}, testfunctionmap = x->x×n) where {D}
#     doubleint = CompDoubleInt(B-> testfunctionmap(B),(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
#     return doubleint
# end

# """
#     volumetrace(Potential)

# Mapping the potential to a volume operator.
# """
# function volumetrace(a::PotentialIntegralOperator{D}, testfunctionmap = x->x) where {D}
#     CompDoubleInt(B->testfunctionmap(B),(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
# end

# trace(pots::LinearCombinationPotentialOperator,interior::Bool;testfunctionmap = x->x) = sum(c*trace(pot,interior;testfunctionmap) for (pot,c) in zip(pots.pots,pots.coeffs))
# ntrace(pots::LinearCombinationPotentialOperator,interior::Bool;testfunctionmap = x->x*n) = sum(c*ntrace(pot,interior;testfunctionmap) for (pot,c) in zip(pots.pots,pots.coeffs))
# rtrace(pots::LinearCombinationPotentialOperator,interior::Bool;testfunctionmap = x->x×n) = sum(c*rtrace(pot,interior;testfunctionmap) for (pot,c) in zip(pots.pots,pots.coeffs))
# ttrace(pots::LinearCombinationPotentialOperator,interior::Bool;testfunctionmap = x->(n× x)×n,testfunction_tangential=false) = sum(c*ttrace(pot,interior;testfunctionmap,testfunction_tangential) for (pot,c) in zip(pots.pots,pots.coeffs))
# pvttrace(pots::LinearCombinationPotentialOperator;testfunctionmap = x->n×(x×n),testfunction_tangential=false) = sum(c*pvttrace(pot;testfunctionmap,testfunction_tangential) for (pot,c) in zip(pots.pots,pots.coeffs))
# pvntrace(pots::LinearCombinationPotentialOperator;testfunctionmap = x->n*x) = sum(c*pvntrace(pot;testfunctionmap) for (pot,c) in zip(pots.pots,pots.coeffs))
# pvtrace(pots::LinearCombinationPotentialOperator;testfunctionmap = x->x) = sum(c*pvtrace(pot;testfunctionmap) for (pot,c) in zip(pots.pots,pots.coeffs))
# pvrtrace(pots::LinearCombinationPotentialOperator;testfunctionmap = x->x×n) = sum(c*pvrtrace(pot;testfunctionmap) for (pot,c) in zip(pots.pots,pots.coeffs))
# volumetrace(pots::LinearCombinationPotentialOperator;testfunctionmap = x->x) = sum(c*volumetrace(pot;testfunctionmap) for (pot,c) in zip(pots.pots,pots.coeffs))

##### from here new trace concept

abstract type AbstractTraceOperator end
struct TraceOperator{D} <: AbstractTraceOperator
    testfunctionmap
    interior::Bool
    principalvalue::Bool
end
TraceOperator{D}(f,interior::Bool) where {D} = TraceOperator{D}(f,interior,false)
function (to::TraceOperator{2})(a::PotentialIntegralOperator{D}) where {D}

    t,f = _trace(a.kernel,to.interior,Val{D}())
    doubleint = CompDoubleInt(B->to.testfunctionmap(B),(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
    to.principalvalue && return doubleint
    t!=0 && (singleint = f*CompSingleInt(B->to.testfunctionmap(B),(x,y)->transpose(x)*y,t,a.op2,a.bfunc))

    if t==0
        return doubleint
    else
        #eventualy a scan to be non zero can be added here.
        return doubleint + singleint
    end

end
function (to::TraceOperator{3})(a::PotentialIntegralOperator{D}) where {D}


    doubleint = CompDoubleInt(B->to.testfunctionmap(B),(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)

    return doubleint


end
function (t::TraceOperator)(a::LinearCombinationPotentialOperator)
    return sum(c*t(a) for (a,c) in zip(a.pots,a.coeffs))
end

const inttrace = TraceOperator{2}(x->x,true)
const exttrace = TraceOperator{2}(x->x,false)
const intntrace = TraceOperator{2}(x->x*n,true)
const extntrace = TraceOperator{2}(x->x*n,false)
const intrtrace = TraceOperator{2}(x->x×n,true)
const extrtrace = TraceOperator{2}(x->x×n,false)
const intttrace = TraceOperator{2}(x->(n× x)×n,true)
const extttrace = TraceOperator{2}(x->(n× x)×n,false)

const pvttrace = TraceOperator{2}(x->n×(x×n),true,true)
const pvntrace = TraceOperator{2}(x->n*x,true,true)
const pvtrace = TraceOperator{2}(x->x,true,true)

"""
    pvrtrace(Potential,interior::Bool)

Compute the principal value of the rotated trace, n× Potential.

This function assumes the normalvector on the mesh to point outwards. 
"""
const pvrtrace = TraceOperator{2}(x->x×n,true,true)

"""
    volumetrace(Potential)

Mapping the potential to a volume operator.
"""
const volumetrace = TraceOperator{3}(x->x,true,true)



## Trace to Direct product space
struct DirectProductTraceOperatorTerm
    traceop 
    hilbertvector
end
struct DirectProductTraceOperator <: AbstractTraceOperator
    terms
    coeffs 
end
Base.getindex(a::AbstractTraceOperator,i::Variational.HilbertVector) = DirectProductTraceOperator([DirectProductTraceOperatorTerm(a,i)],[1.0])
+(a::DirectProductTraceOperator,b::DirectProductTraceOperator) = DirectProductTraceOperator([a.terms;b.terms],[a.coeffs;b.coeffs])
*(a::Number,b::DirectProductTraceOperator) = DirectProductTraceOperator(b.terms,a*b.coeffs)
-(a::DirectProductTraceOperator) = DirectProductTraceOperator(a.terms,-a.coeffs)

function (t::DirectProductTraceOperator)(a::DirectProductPotentialOperator)
    return sum(c1*c2*traceterm.traceop(potterm.pot)[traceterm.hilbertvector,potterm.hilbertvector] 
    for (traceterm,c1) in zip(t.terms,t.coeffs) for (potterm,c2) in zip(a.terms,a.coeffs))
end
