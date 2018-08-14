"""
Solves a variational equation by simply creating the full system matrix
and calling a traditional lu decomposition.
"""
function solve(eq)

    time_domain = isa(first(eq.trial_space_dict).second, BEAST.SpaceTimeBasis)
    time_domain |= isa(first(eq.trial_space_dict).second, BEAST.StagedTimeStep)
    if time_domain
        return td_solve(eq)
    end

    test_space_dict  = eq.test_space_dict
    trial_space_dict = eq.trial_space_dict

    lhs = eq.equation.lhs
    rhs = eq.equation.rhs

    b = assemble(rhs, test_space_dict)
    Z = assemble(lhs, test_space_dict, trial_space_dict)

    u = Z \ b

    return u
end


function td_solve(eq)

    @warn "very limited support for automated solution of TD equations...."
    op = eq.equation.lhs.terms[1].kernel
    fn = eq.equation.rhs.terms[1].functional

    V = eq.trial_space_dict[1]
    W = eq.test_space_dict[1]

    A = assemble(op, W, V)
    S = inv(A[:,:,1])
    b = assemble(fn, W)

    nt = numfunctions(temporalbasis(V))
    marchonintime(S, A, b, nt)
end
