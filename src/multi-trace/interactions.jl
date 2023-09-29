

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
    # if Ω2!==Ω3
    #     sign = -1
    # elseif Ω2===Ω3
    #     sign = 1
    # end
    @warn "correct sign for inside is 1?"
    sign = 1
    
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
struct VectorStrat <: BEAST.NumericalStrategy end

function convert_inside_to_outside_basis(child,parent,::VectorStrat)
    a = [-1 0 0 0
      0 1 0 0
      0 0 -parent.data.μ/child.data.μ 0
      0 0 0 -child.data.ϵ/parent.data.ϵ]

      return a
end




function (int::Interaction{<: Domain{HomogeneousDomain},<: Domain{HomogeneousDomain},<: RootDomain})(::VectorStrat)
    k = sqrt(int.embedvol.data.ϵ*int.embedvol.data.μ)*int.embedvol.data.ω
    green = HHH.green(wavenumber=k)
    gradgreen = HHH.gradgreen(wavenumber=k)
    b = basisfunction()
    @warn "check if extra - in front of a is correct, describtion of As from paper asumes n outward so inward in outer domain?"
    a = -[n×(gradgreen×nothing)          n×(green(n*b))        -(n×green)                n×gradgreen
        BEAST.ZeroOperator()                   -(gradgreen⋅nothing)(n*b)     (gradgreen⋅nothing)               -(-k^2*green)
        -(n×(gradgreen(∇⋅b)))-k^2*(n×green)    -(n×((gradgreen×nothing)(n*b)))    n×(gradgreen×nothing)   BEAST.ZeroOperator()
        -(n⋅(gradgreen×nothing))              -(n⋅green(n*b))         n⋅green -(n⋅gradgreen)]
    id = [Identity() ZeroOperator() ZeroOperator() ZeroOperator()
        ZeroOperator() Identity() ZeroOperator() ZeroOperator()
        ZeroOperator() ZeroOperator() Identity() ZeroOperator()
        ZeroOperator() ZeroOperator() ZeroOperator() Identity()]
    a = id - a
    
    if (int.testvol.id,int.trialvol.id) in keys(int.config.touching) 
        println("cauchy limit taken")
        a = BEAST.cauchylimit.(a;Ω1=int.testvol,Ω2=int.trialvol,Ω3=int.embedvol)
    end
    a = BEAST.normalorient.(a;Ω1=int.testvol,Ω2=int.trialvol,Ω3=int.embedvol)
    
    return a

end
function (int::Interaction{<: Domain{HomogeneousDomain},<: Domain{HomogeneousDomain},<: SubDomain})(::VectorStrat)
    k = sqrt(int.embedvol.data.ϵ*int.embedvol.data.μ)*int.embedvol.data.ω #foute lijn, puur voor test!!!
    green = HHH.green(wavenumber=k)
    gradgreen = HHH.gradgreen(wavenumber=k)
    b = basisfunction()

    a = -[n×(gradgreen×nothing)          n×(green(n*b))        -(n×green)                n×gradgreen
        BEAST.ZeroOperator()                   -(gradgreen⋅nothing)(n*b)     (gradgreen⋅nothing)               -(-k^2*green)
        -(n×(gradgreen(∇⋅b)))-k^2*(n×green)    -(n×((gradgreen×nothing)(n*b)))    n×(gradgreen×nothing)   BEAST.ZeroOperator()
        -(n⋅(gradgreen×nothing))              -(n⋅green(n*b))         n⋅green -(n⋅gradgreen)]
        id = [Identity() ZeroOperator() ZeroOperator() ZeroOperator()
        ZeroOperator() Identity() ZeroOperator() ZeroOperator()
        ZeroOperator() ZeroOperator() Identity() ZeroOperator()
        ZeroOperator() ZeroOperator() ZeroOperator() Identity()]
    a = id - a
    
    if (int.testvol.id,int.trialvol.id) in keys(int.config.touching) 
        println("cauchy limit taken")
        a = BEAST.cauchylimit.(a;Ω1=int.testvol,Ω2=int.trialvol,Ω3=int.embedvol)
    end
   a = BEAST.normalorient.(a;Ω1=int.testvol,Ω2=int.trialvol,Ω3=int.embedvol)
    return a

end