import Base:*, +, -

mutable struct DiscreteEq 
    lhs
    rhs
    trialspace
end

function +(a::VGL,b::VGL)
    @assert a.trialspace == b.trialspace
    return VGL(a.lhs+b.lhs,a.rhs+b.rhs,a.trialspace)
end

*(a::Union{LinearMap,Number},b::VGL) = VGL(a*b.lhs,a*b.rhs,b.trialspace)
-(a::VGL,b::VGL) = a + (-1)*b

function generate_discrete_equations(w::World,strat;kwargs...)
    # create hilbertspaces


end

### to define for your specific problem
function hilbertspace(obj,strat)
    return true
end

function getequations(strat,dict)# builds the equations in vp case there are 2 homogeneous equations and 2 pec cfie equations
    
    return []
end
##### vanaf hier voor vp

function getequations(strat::VPPMCHWT,d::Dict)
    homind = d[HOM()]
    mfie = d[MFIE()]
    efie = d[EFIE()]
    cfie = d[CFIE()]
    freespace = d[FREESPACE()]
    eq1 = []
end
function testspace(obj,strat) end
function trialspace(obj,strat) end
function getequations(strat::VP,dict) end