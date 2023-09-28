

"""
Write for each type of interaction this type of function.
"""
function (int::Interaction{<:Domain{BackgroundDomain},})() 
    i = [] #interaction matrix
    normalorient.(i,Ω1,Ω2,Ω3) # example 

    return nothing
end

function alpha(Ω1,Ω2)


end
"""
cauchylimit needs to be applied before the normalorient function
"""

function cauchylimit(operator::AbstractOperator;Ω1,Ω2,Ω3)
#check first if touching is non empty
@assert is_child_of(Ω1,Ω3)||Ω1===Ω3
@assert is_child_of(Ω2,Ω3)||Ω2===Ω3
    if Ω2!==Ω3
        sign = -1
    elseif Ω2===Ω3
        sign = 1
    end
    
    trace(operator,sign)

end 

function trace(op::AbstractOperator,sign)
    @warn "general abstract opterator trace function called returning pv of operator!: "*string(typeof(op))
    return op
end
function normalorient(op::AbstractOperator,signtest,signtrial)
    @warn "normalorient not implemented for: "*string(typeof(op))
    return op
end
trace(op::ZeroOperator,s) = op

function trace(op::LinearCombinationOfOperators,sign)
    result = ZeroOperator()
    for (c,o) in zip(op.coeffs,op.ops)
        result += c*trace(o,sign)
    end
    return result
end


function normalorient(op::LinearCombinationOfOperators,signtest,signtrial)
    result = ZeroOperator()
    for (c,o) in zip(op.coeffs,op.ops)
        result += c*normalorient(o,signtest,signtrial)
    end
    return result
end

function normalorient(operator::AbstractOperator;Ω1,Ω2,Ω3) 
    if Ω1===Ω3
        sign_test_normal = 1
    else
        sign_test_normal = -1
    end
    if Ω2===Ω3
        sign_trial_normal = 1
    else
        sign_trial_normal = -1
    end
    normalorient(operator,sign_test_normal,sign_trial_normal)
end

normalorient(op::ZeroOperator,a,b) = op

###### Interactions