using LinearAlgebra

# function convert_to_dense(A::LinearMaps.LinearCombination)

#     # @info "convert matrix to dense..."
#     # T = eltype(A)
#     # # M = zeros(T, size(A))
#     # # M = PseudoBlockMatrix{T}(undef, BlockArrays.blocksizes(A)...)
#     # # fill!(M,0)
#     # M = zeros(eltype(A), axes(A))
#     # @show typeof(M)
#     # for map in A.maps
#     #     mul!(M, 1, map, 1, 1)
#     # end

#     # @info "matrix converted to dense"
#     # return M
#     return Matrix(A)
# end

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

    X = _spacedict_to_directproductspace(eq.test_space_dict)
    Y = _spacedict_to_directproductspace(eq.trial_space_dict)

    b = assemble(rhs, X)
    Z = assemble(lhs, X, Y)
    # b = assemble(rhs, test_space_dict)
    # Z = assemble(lhs, test_space_dict, trial_space_dict)
    # M = convert_to_dense(Z)
    # println("Sysmatrix converted to dense.")

    print("Converting system to Matrix...")
    M = Matrix(Z)
    println("done.")

    print("LU solution of the linear system...")
    u = M \ Vector(b)
    println("done.")

    # return u
    # return PseudoBlockVector(u, blocksizes(Z,2))
    ax = nestedrange(Y, 1, numfunctions)
    return PseudoBlockVector(u, (ax,))
end


function td_solve(eq)

    V = eq.trial_space_dict[1]

    # bilform = eq.equation.lhs

    A = td_assemble(eq.equation.lhs, eq.test_space_dict, eq.trial_space_dict)
    T = eltype(A)
    S = zeros(T, size(A)[1:2])
    ConvolutionOperators.timeslice!(S, A, 1)

    # S = timeslice(A,1)
    # iS = inv(Array(S))
    iS = inv(S)

    b = td_assemble(eq.equation.rhs, eq.test_space_dict)

    nt = numfunctions(temporalbasis(V))
    marchonintime(iS, A, b, nt)
end

# function timeslice(A::BlockArray, k)

#     I = blocklengths(axes(A,1))
#     J = blocklengths(axes(A,2))

#     T = eltype(eltype(A))
#     S = PseudoBlockArray{T}(undef, I, J)
#     fill!(S,0)

#     for i in 1:blocksize(A,1)
#         for j in 1:blocksize(A,2)
#             isassigned(A.blocks, i, j) || continue
#             A[Block(i,j)] isa Zeros && continue
#             A[Block(i,j)] isa Fill && continue
#             Aij = A[Block(i,j)]
#             Ak = AbstractArray(A[Block(i,j)].convop)[:,:,k]
#             S[Block(i,j)] = Ak
#         end
#     end

#     return S
# end
