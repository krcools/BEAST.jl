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

    # @warn "very limited support for automated solution of TD equations...."
    # op = eq.equation.lhs.terms[1].kernel
    # fn = eq.equation.rhs.terms[1].functional

    V = eq.trial_space_dict[1]
    # W = eq.test_space_dict[1]

    bilform = eq.equation.lhs
    # I = [numfunctions(spatialbasis(eq.test_space_dict[i])) for i in 1:length(bilform.test_space)]
    # J = [numfunctions(spatialbasis(eq.trial_space_dict[i])) for i in 1:length(bilform.trial_space)]

    A = td_assemble(eq.equation.lhs, eq.test_space_dict, eq.trial_space_dict)

    # T = eltype(eltype(A))
    # S = PseudoBlockArray{T}(undef, I, J)
    # fill!(S,0)
    #
    # for i in 1:nblocks(A,1)
    #     for j in 1:nblocks(A,2)
    #         isassigned(A.blocks, i, j) || continue
    #         S[Block(i,j)] = A[Block(i,j)].banded[:,:,1]
    #     end
    # end

    S = timeslice(A,1)
    iS = inv(Array(S))

    b = td_assemble(eq.equation.rhs, eq.test_space_dict)

    nt = numfunctions(temporalbasis(V))
    # @show @which marchonintime(AiS, A, b, nt)
    marchonintime(iS, A, b, nt)
end

function timeslice(A::BlockArray, k)

    I = [blocksize(A, (1,i))[1] for i in 1:nblocks(A,1)]
    J = [blocksize(A, (j,1))[2] for j in 1:nblocks(A,2)]

    T = eltype(eltype(A))
    S = PseudoBlockArray{T}(undef, I, J)
    fill!(S,0)

    for i in 1:nblocks(A,1)
        for j in 1:nblocks(A,2)
            isassigned(A.blocks, i, j) || continue
            S[Block(i,j)] = A[Block(i,j)].banded[:,:,k]
        end
    end

    return S
end
