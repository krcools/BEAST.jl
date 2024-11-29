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

"""
    ttrace(Potential,interior::Bool;testfunction_tangential=false)

Compute the tangential trace, -n×n× Potential, of a potential operator mapping it to a boundary operator.

This function assumes the normalvector on the mesh to point outwards. 
Global multi-trace structures are supported when this function is called.
A small acceleration can be gained by setting testfunction_tangential to true.
"""
function ttrace(a::PotentialIntegralOperator{D},interior::Bool;testfunction_tangential=false) where {D}
    t,f = _trace(a.kernel,interior,Val{D}())
    if testfunction_tangential
        doubleint = CompDoubleInt(B->B,(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
        t!=0 && (singleint = f*CompSingleInt(B->B,(x,y)->transpose(x)*y,t,a.op2,a.bfunc))
    else

        doubleint = -1*CompDoubleInt(B-> n×(n×B),(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
        t!=0 && (singleint = (-f)*CompSingleInt(B -> n×(n×B),(x,y)->transpose(x)*y,t,a.op2,a.bfunc))
    end
    if t==0
     
        return doubleint
    else
       
        #eventualy a scan to be non zero can be added here.
        return doubleint + singleint
    end
end
"""
    ntrace(Potential,interior::Bool)

Compute the normal trace, n⋅ Potential, of a potential operator mapping it to a boundary operator.

This function assumes the normalvector on the mesh to point outwards. 
Global multi-trace structures are supported when this function is called.
"""
function ntrace(a::PotentialIntegralOperator{D},interior::Bool) where {D}
    t,f = _trace(a.kernel,interior,Val{D}())
    doubleint = CompDoubleInt(B-> B*n,(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
    t!=0 && (singleint = f*CompSingleInt(B-> B*n,(x,y)->transpose(x)*y,t,a.op2,a.bfunc))
    
    if t==0
        return doubleint
    else
        #eventualy a scan to be non zero can be added here.
        return doubleint + singleint
    end
end
"""
    trace(Potential,interior::Bool)

Compute the trace of a potential operator mapping it to a boundary operator.

This function assumes the normalvector on the mesh to point outwards. 
Global multi-trace structures are supported when this function is called.
"""
function trace(a::PotentialIntegralOperator{D},interior::Bool) where {D}
    t,f = _trace(a.kernel,interior,Val{D}())
    doubleint = CompDoubleInt(B->B,*,a.kernel,a.op2,a.bfunc)
    t!=0 && (singleint = f*CompSingleInt(B->B,*,t,a.op2,a.bfunc))
    
    if t==0
        return doubleint
    else
        #eventualy a scan to be non zero can be added here.
        return doubleint + singleint
    end
end
"""
    rtrace(Potential,interior::Bool)

Compute the rotated trace, n × Potential, of a potential operator mapping it to a boundary operator.

This function assumes the normalvector on the mesh to point outwards. 
Global multi-trace structures are supported when this function is called.
"""
function rtrace(a::PotentialIntegralOperator{D},interior::Bool) where {D}
    t,f = _trace(a.kernel,interior,Val{D}())
    doubleint = CompDoubleInt(B-> B × n,(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
    t!=0 && (singleint = f*CompSingleInt(B-> B × n,(x,y)->transpose(x)*y,t,a.op2,a.bfunc))
    
    if t==0
        return doubleint
    else
        #eventualy a scan to be non zero can be added here.
        return doubleint + singleint
    end
end
"""
    pvttrace(Potential,interior::Bool;testfunction_tangential=false)

Compute the principal value of the tangential trace, -n×n× Potential.

A small acceleration can be gained by setting testfunction_tangential to true.
"""
# the principal value of the trace.
function pvttrace(a::PotentialIntegralOperator{D};testfunction_tangential=false) where {D}
    if testfunction_tangential
        doubleint = CompDoubleInt(B->B,(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
    else
        doubleint = -1*CompDoubleInt(B-> n×(n×B),(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
    end
    return doubleint
end

"""
    pvntrace(Potential,interior::Bool)

Compute the principal value of the normal trace, n⋅ Potential.

This function assumes the normalvector on the mesh to point outwards. 
"""
function pvntrace(a::PotentialIntegralOperator{D}) where {D}

    doubleint = CompDoubleInt(B-> B*n,(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)

    return doubleint
end

"""
    pvtrace(Potential,interior::Bool)

Compute the principal value of the  trace.
"""
function pvtrace(a::PotentialIntegralOperator{D}) where {D}
    doubleint = CompDoubleInt(B-> B,*,a.kernel,a.op2,a.bfunc)

    return doubleint
end

"""
    pvrtrace(Potential,interior::Bool)

Compute the principal value of the rotated trace, n× Potential.

This function assumes the normalvector on the mesh to point outwards. 
"""
function pvrtrace(a::PotentialIntegralOperator{D}) where {D}
    doubleint = CompDoubleInt(B-> B × n,(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
    return doubleint
end

"""
    volumetrace(Potential)

Mapping the potential to a volume operator.
"""
function volumetrace(a::PotentialIntegralOperator{D}) where {D}
    CompDoubleInt(B->B,(x,y)->transpose(x)*y,a.kernel,a.op2,a.bfunc)
end





