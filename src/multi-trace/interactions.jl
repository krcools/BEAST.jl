
"""
A is the type of the test structure, B is the type of the trial structure, C the type of ths structure embedding
"""
struct Interaction{A,B,C}# wordt gecrieerd bij de itteratie over alle volumes.
    config
    testvol::A
    trialvol::B
    embedvol::C
end
"""
Write for each type of interaction this type of function.
"""
function (int::Interaction)() 
    i = [] #interaction matrix
    normalorient.(i,Ω1,Ω2,Ω3) # example 

    return nothing
end

function alpha(Ω1,Ω2)


end
"""
cauchylimit needs to be applied before the normalorient function
"""

function cauchylimit(operator::AbstractOperator,Ω1,Ω2,Ω3)
#check first if touching is non empty
@assert is_child_of(Ω1,Ω3)||Ω1===Ω3
@assert is_child_of(Ω2,Ω3)||Ω2===Ω3
    if !Ω2===Ω3
        sign = -1
    elseif Ω2===Ω3
        sign = 1
    end
    
    trace(operator,sign)

end 

function trace(op::AbstractOperator,sign)
    @warning "general abstract opterator trace function called returning pv of operator!: "*string(typeof(op))
    return op
end
function normalorient(op::AbstractOperator,signtest,signtrial)
    @warning "normalorient not implemented for: "*string(typeof(op))
    return op
end

function trace(op::LinearCombinationOfOperators,sign)
    result = ZeroOperator()
    for (c,o) in zip(op.coeffs,op.ops)
        result += c*trace(o,sign)
    end
    return result
end


function normalorient(operator::LinearCombinationOfOperators,signtest,signtrial)
    result = ZeroOperator()
    for (c,o) in zip(op.coeffs,op.ops)
        result += c*normalorient(o,signtest,signtrial)
    end
    return result
end

function normalorient(operator::AbstractOperator,Ω1,Ω2,Ω3) 
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


###### Interactions