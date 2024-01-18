using LinearAlgebra

lusolve(eq) = solve(eq)

"""
Solves a variational equation by simply creating the full system matrix
and calling a traditional lu decomposition.
"""
function solve(eq)

    time_domain = isa(first(eq.trial_space_dict).second, BEAST.SpaceTimeBasis)
    time_domain_cq = isa(first(eq.trial_space_dict).second, BEAST.StagedTimeStep)
    if time_domain
        return td_solve(eq)
    end

    if time_domain_cq
        return td_solve_cq(eq)
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
    return PseudoBlockVector(u, (ax,))
end

function td_solve_cq(eq)

    @warn("very limited sypport for automated solution of TD equations....")
    op = eq.equation.lhs.terms[1].kernel
    fn = eq.equation.rhs.terms[1].functional

    V = eq.trial_space_dict[1]
    W = eq.test_space_dict[1]

    A = assemble(op, W, V)
    S = inv(A[:,:,1])
    b = assemble(fn, W)

    nt = numfunctions(temporalbasis(V))
    marchonintime_cq(S, A, b, nt)
end
