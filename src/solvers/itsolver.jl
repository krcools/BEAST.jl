import IterativeSolvers




struct GMRESSolver{L,T} <: LinearMap{T}
    linear_operator::L
    maxiter::Int
    restart::Int
    tol::T
end


Base.axes(A::GMRESSolver) = reverse(axes(A.linear_operator))


function GMRESSolver(op; maxiter=0, restart=0, tol=sqrt(eps(real(eltype(op)))))

    m, n = size(op)
    @assert m == n

    maxiter == 0 && (maxiter = div(n, 5))
    restart == 0 && (restart = n)

    GMRESSolver(op, maxiter, restart, tol)
end


operator(solver::GMRESSolver) = solver.linear_operator


# function solve(solver::GMRESSolver, b)
#     op = operator(solver)
#     x, ch = IterativeSolvers.gmres(op, b, log=true,  maxiter=solver.maxiter,
#         restart=solver.restart, reltol=solver.tol, verbose=true)
#     return x, ch
# end


function solve!(x, solver::GMRESSolver, b)
    op = operator(solver)
    x, ch = IterativeSolvers.gmres!(x, op, b, log=true,  maxiter=solver.maxiter,
        restart=solver.restart, reltol=solver.tol, verbose=true)
    return x, ch
end


function Base.:*(A::GMRESSolver, b::AbstractVector)

    T = promote_type(eltype(A), eltype(b))
    y = PseudoBlockVector{T}(undef, BlockArrays.blocksizes(A)[1])

    mul!(y, A, b)

    # x, ch = solve(solver, b)
    # println("Number of iterations: ", ch.iters)
    # ch.isconverged || error("Iterative solver did not converge.")
    # return x
end

Base.size(solver::GMRESSolver) = reverse(size(solver.linear_operator))

function LinearAlgebra.mul!(y::AbstractVecOrMat, solver::GMRESSolver, x::AbstractVector)
    fill!(y,0)
    y, ch = solve!(y, solver, x)
    println("Number of iterations: ", ch.iters)
    ch.isconverged || error("Iterative solver did not converge.")
    return y
end


function gmres(eq::DiscreteEquation; maxiter=0, restart=0, tol=0)

    test_space_dict  = eq.test_space_dict
    trial_space_dict = eq.trial_space_dict

    lhs = eq.equation.lhs
    rhs = eq.equation.rhs

    b = assemble(rhs, test_space_dict)
    Z = assemble(lhs, test_space_dict, trial_space_dict)

    block_sizes = zeros(Int, length(trial_space_dict))
    for (p,x) in eq.trial_space_dict
        block_sizes[p] = numfunctions(x)
    end

    T = promote_type(eltype(Z), eltype(b))
    x = PseudoBlockVector{T}(undef, block_sizes)
    fill!(x, 0)

    if tol == 0
        invZ = GMRESSolver(Z, maxiter=maxiter, restart=restart)
    else
        invZ = GMRESSolver(Z, maxiter=maxiter, restart=restart, tol=tol)
    end
    # x = invZ * b
    mul!(x, invZ, b)
end
