using LinearAlgebra

lusolve(eq) = solve(eq)

"""
Solves a variational equation by simply creating the full system matrix
and calling a traditional lu decomposition.
"""
function solve(eq;Zz=false)

    time_domain = isa(first(eq.trial_space_dict).second, BEAST.SpaceTimeBasis)
    time_domain |= isa(first(eq.trial_space_dict).second, BEAST.StagedTimeStep)
    if time_domain
        return td_solve(eq)
    end

    test_space_dict  = eq.test_space_dict
    trial_space_dict = eq.trial_space_dict

    lhs = eq.equation.lhs
    rhs = eq.equation.rhs

    X = _spacedict_to_directproductspace(eq.test_space_dict)
    Y = _spacedict_to_directproductspace(eq.trial_space_dict)

    b = assemble(rhs, X)
    Z = assemble(lhs, X, Y)

    print("Converting system to Matrix...")
    M = Matrix(Z)
    println("done.")

    print("LU solution of the linear system...")
    u = M \ Vector(b)
    println("done.")

    ax = nestedrange(Y, 1, numfunctions)
    Zz && (return PseudoBlockVector(u, (ax,)), Z)
    return PseudoBlockVector(u, (ax,))
end


function solve(b,Z,Y;Zz=true)

    
    print("Converting system to Matrix...")
    M = Matrix(Z)
    println("done.")

    print("LU solution of the linear system...")
    u = M \ Vector(b)
    println("done.")

    ax = nestedrange(Y, 1, numfunctions)
    Zz && (return PseudoBlockVector(u, (ax,)), Z)
    return PseudoBlockVector(u, (ax,))
end
