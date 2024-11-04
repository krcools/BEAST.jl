
#### trace vector 

# orientation is used when test and trial chart have same displacement, in this case the test-cart is displaced in the orientation*normal(test-chart) direction.
struct DisplacementVector{T} <: Kernel{T} 
    orientation
end


function (u::DisplacementVector)(x,disp)
    if disp == 0
        return sign(u.orientation)*normal(x)
    else
        return sign(disp)*normal(x)
    end
end


#### trace of a potential:
_trace(::Kernel,n;type=Float64,D) = 0, 0
function _trace(::HH3DGradGreen,n,D;type=Float64)
    if D==2
        return DisplacementVector{type}(n) , 0.5
    end
    return 0, 0
end

"""
sign indicates if trace is taken allong the normal (+1) are opposite to the normal (-1)
"""
function ttrace(a::PotentialIntegralOperator{D},sign;testfunction_tangential=false) where {D}
    t,f = _trace(a.kernel,sign,D)
    if testfunction_tangential
        doubleint = CompDoubleInt(B->B,(x,y)->transpose(x)*y,a.kernel,a.pairing2,a.bfunc)
        t!=0 && singleint = f*CompSingleInt(B->B,(x,y)->transpose(x)*y,t,a.pairing2,a.bfunc)
    else

        doubleint = -1*CompDoubleInt(B-> n×(n×B),(x,y)->transpose(x)*y,a.kernel,a.pairing2,a.bfunc)
        t!=0 && singleint = (-f)*CompSingleInt(B -> n×(n×B),(x,y)->transpose(x)*y,t,a.pairing2,a.bfunc)
    end
    if t==0
        return doubleint
    else
        #eventualy a scan to be non zero can be added here.
        return doubleint + singleint
    end
end

function ntrace(a::PotentialIntegralOperator{D},sign) where {D}
    t,f = _trace(a.kernel,sign,D)
    doubleint = CompDoubleInt(B-> B*n,(x,y)->transpose(x)*y,a.kernel,a.pairing2,a.bfunc)
    t!=0 && singleint = f*CompSingleInt(B-> B*n,(x,y)->transpose(x)*y,t,a.pairing2,a.bfunc)
    
    if t==0
        return doubleint
    else
        #eventualy a scan to be non zero can be added here.
        return doubleint + singleint
    end
end

function trace(a::PotentialIntegralOperator{D},sign) where {D}
    t,f = _trace(a.kernel,sign,D)
    doubleint = CompDoubleInt(B->B,*,a.kernel,a.pairing2,a.bfunc)
    t!=0 && singleint = f*CompSingleInt(B->B,*,t,a.pairing2,a.bfunc)
    
    if t==0
        return doubleint
    else
        #eventualy a scan to be non zero can be added here.
        return doubleint + singleint
    end
end

function strace(a::PotentialIntegralOperator{D},sign) where {D}
    t,f = _trace(a.kernel,sign,D)
    doubleint = CompDoubleInt(B-> B × n,(x,y)->transpose(x)*y,a.kernel,a.pairing2,a.bfunc)
    t!=0 && singleint = f*CompSingleInt(B-> B × n,(x,y)->transpose(x)*y,t,a.pairing2,a.bfunc)
    
    if t==0
        return doubleint
    else
        #eventualy a scan to be non zero can be added here.
        return doubleint + singleint
    end
end

# the principal value of the trace.
function pvttrace(a::PotentialIntegralOperator{D},sign;testfunction_tangential=false) where {D}
    if testfunction_tangential
        doubleint = CompDoubleInt(B->B,(x,y)->transpose(x)*y,a.kernel,a.pairing2,a.bfunc)
    else
        doubleint = -1*CompDoubleInt(B-> n×(n×B),(x,y)->transpose(x)*y,a.kernel,a.pairing2,a.bfunc)
    end
    return doubleint
end

function pvntrace(a::PotentialIntegralOperator{D},sign) where {D}

    doubleint = CompDoubleInt(B-> B*n,(x,y)->transpose(x)*y,a.kernel,a.pairing2,a.bfunc)

    return doubleint
end

function pvtrace(a::PotentialIntegralOperator{D},sign) where {D}
    doubleint = CompDoubleInt(B-> B,*,a.kernel,a.pairing2,a.bfunc)

    return doubleint
end

function pvstrace(a::PotentialIntegralOperator{D},sign) where {D}
    doubleint = CompDoubleInt(B-> B × n,(x,y)->transpose(x)*y,a.kernel,a.pairing2,a.bfunc)
    return doubleint
end

function volumetrace(a::PotentialIntegralOperator{D}) where {D}
    CompDoubleInt(B->B,(x,y)->transpose(x)*y,a.kernel,a.pairing2,a.bfunc)
end





