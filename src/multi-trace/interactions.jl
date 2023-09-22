
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

function cauchylimit(operator,Ω1,Ω2,Ω3)
#check first if touching is non empty
@assert is_child_of(Ω1,Ω3)||Ω1===Ω3
@assert is_child_of(Ω2,Ω3)||Ω2===Ω3
    if !Ω2===Ω3
        sign = -1
    elseif Ω2===Ω3
        sign = 1
    end
    
    trace(operator,sign)

end #TODO define these


function normalorient(operator,Ω1,Ω2,Ω3) end


###### Interactions