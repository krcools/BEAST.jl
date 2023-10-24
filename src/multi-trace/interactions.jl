

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

function convert_outside_to_inside_basis(child,parent,::VectorStrat)
    p = physicalconstants(parent.data)
    c = physicalconstants(child.data)
    a = [-1 0 0 0
      0 1 0 0
      0 0 -c.μ/p.μ 0
      0 0 0 -p.ϵ/c.ϵ]

      return a
end




# function (int::Interaction{<: Domain{HomogeneousDomain},<: Domain{HomogeneousDomain},<: RootDomain})(::VectorStrat)
#     p = physicalconstants(int.embedvol.data)
#     k = sqrt(p.ϵ*p.μ)*p.ω
#     green = HHH.green(wavenumber=k)
#     gradgreen = HHH.gradgreen(wavenumber=k)
#     b = basisfunction()
#     @warn "check if extra - in front of a is correct, describtion of As from paper asumes n outward so inward in outer domain?"
#     a = -[n×(gradgreen×nothing)          n×(green(n*b))        -(n×green)                n×gradgreen
#         BEAST.ZeroOperator()                   -(gradgreen⋅nothing)(n*b)     (gradgreen⋅nothing)               -(-k^2*green)
#         -(n×(gradgreen(∇⋅b)))-k^2*(n×green)    -(n×((gradgreen×nothing)(n*b)))    n×(gradgreen×nothing)   BEAST.ZeroOperator()
#         -(n⋅(gradgreen×nothing))              -(n⋅green(n*b))         n⋅green -(n⋅gradgreen)]
#     id = [Identity() ZeroOperator() ZeroOperator() ZeroOperator()
#         ZeroOperator() Identity() ZeroOperator() ZeroOperator()
#         ZeroOperator() ZeroOperator() Identity() ZeroOperator()
#         ZeroOperator() ZeroOperator() ZeroOperator() Identity()]
#     #a = id - a
    
#     if (int.testvol.id,int.trialvol.id) in keys(int.config.touching) 
#         println("cauchy limit taken")
#         a = BEAST.cauchylimit.(a;Ω1=int.testvol,Ω2=int.trialvol,Ω3=int.embedvol)
#     end
#     a = BEAST.normalorient.(a;Ω1=int.testvol,Ω2=int.trialvol,Ω3=int.embedvol)
    
#     return a

# end
# function (int::Interaction{<: Domain{HomogeneousDomain},<: Domain{HomogeneousDomain},<: SubDomain})(::VectorStrat)
#     p = physicalconstants(int.embedvol.data)
#     k = sqrt(p.ϵ*p.μ)*p.ω
#     green = HHH.green(wavenumber=k)
#     gradgreen = HHH.gradgreen(wavenumber=k)
#     b = basisfunction()

#     a = -[n×(gradgreen×nothing)          n×(green(n*b))        -(n×green)                n×gradgreen
#         BEAST.ZeroOperator()                   -(gradgreen⋅nothing)(n*b)     (gradgreen⋅nothing)               -(-k^2*green)
#         -(n×(gradgreen(∇⋅b)))-k^2*(n×green)    -(n×((gradgreen×nothing)(n*b)))    n×(gradgreen×nothing)   BEAST.ZeroOperator()
#         -(n⋅(gradgreen×nothing))              -(n⋅green(n*b))         n⋅green -(n⋅gradgreen)]
#         id = [Identity() ZeroOperator() ZeroOperator() ZeroOperator()
#         ZeroOperator() Identity() ZeroOperator() ZeroOperator()
#         ZeroOperator() ZeroOperator() Identity() ZeroOperator()
#         ZeroOperator() ZeroOperator() ZeroOperator() Identity()]
#     #a = id - a
    
#     if (int.testvol.id,int.trialvol.id) in keys(int.config.touching) 
#         println("cauchy limit taken")
#         a = BEAST.cauchylimit.(a;Ω1=int.testvol,Ω2=int.trialvol,Ω3=int.embedvol)
#     end
#    a = BEAST.normalorient.(a;Ω1=int.testvol,Ω2=int.trialvol,Ω3=int.embedvol)
#     return a

# end


function (int::Interaction{<: Domain{HomogeneousDomain},<: Domain{HomogeneousDomain},<:Union{RootDomain,SubDomain}})(::VectorStrat)
    p = physicalconstants(int.embedvol.data)
    k = sqrt(p.ϵ*p.μ)*p.ω
    G = BEAST.greenhh3d(wavenumber=k)
    ∇G = BEAST.∇(G)
    Ω1=int.testvol
    Ω2=int.trialvol
    Ω3=int.embedvol
    ∇Gx =  BEAST.build_potential(∇G×B,Ω2.data.trialbasises[1].geo)
    Gn = BEAST.build_potential(G*n*B,Ω2.data.trialbasises[1].geo)
    Gnx = BEAST.build_potential(G*n × B,Ω2.data.trialbasises[1].geo)
    ∇G = BEAST.build_potential(∇G*B,Ω2.data.trialbasises[1].geo)
    ∇Gdotn = BEAST.build_potential(∇G⋅n*B,Ω2.data.trialbasises[1].geo)
    ∇Gdot = BEAST.build_potential(∇G⋅B,Ω2.data.trialbasises[1].geo)
    
    Gr = BEAST.build_potential(G*B,Ω2.data.trialbasises[1].geo)
    ∇G∇B = BEAST.build_potential(∇G*∇(B),Ω2.data.trialbasises[1].geo)
    ∇Gxn = BEAST.build_potential(∇G×n*B,Ω2.data.trialbasises[1].geo)
    testsurf = Ω1.data.testbasises[1].geo
    if Ω1==Ω2
    a = -[γₛ(∇Gx,testsurf)         γₛ(Gn,testsurf)      -γₛ(Gnx,testsurf)           γₛ(∇G,testsurf)
        BEAST.ZeroOperator()            -τ(∇Gdotn,testsurf)    τ(∇Gdot,testsurf)     k^2*τ(Gr,testsurf)
        -γₛ(∇G∇B,testsurf)-k^2*γₛ(Gr,testsurf)   -γₛ(∇Gxn,testsurf)  γₛ(∇Gx,testsurf)  BEAST.ZeroOperator()
         -γₙ(∇Gx,testsurf)         -γₙ(Gn,testsurf)     γₙ(Gr,testsurf)  -γₙ(∇G,testsurf)]
    else
        a = -[γₛᶜ(∇Gx,testsurf)         γₛᶜ(Gn,testsurf)      -γₛᶜ(Gnx,testsurf)           γₛᶜ(∇G,testsurf)
        BEAST.ZeroOperator()            -τᶜ(∇Gdotn,testsurf)    τᶜ(∇Gdot,testsurf)     k^2*τᶜ(Gr,testsurf)
        -γₛᶜ(∇G∇B,testsurf)-k^2*γₛᶜ(Gr,testsurf)   -γₛᶜ(∇Gxn,testsurf)  γₛᶜ(∇Gx,testsurf)  BEAST.ZeroOperator()
         -γₙᶜ(∇Gx,testsurf)         -γₙᶜ(Gn,testsurf)     γₙᶜ(Gr,testsurf)  -γₙᶜ(∇G,testsurf)]
    end


    a = BEAST.normalorient.(a;Ω1=int.testvol,Ω2=int.trialvol,Ω3=int.embedvol)
    
    return a

end