import IterativeSolvers




struct GMRESSolver{L,R}
    linear_operator::L
    maxiter::Int
    restart::Int
    tol::R
end


function GMRESSolver(op; maxiter=0, restart=0, tol=sqrt(eps(real(eltype(op)))))

    m, n = size(op)
    @assert m == n

    maxiter == 0 && (maxiter = div(n, 5))
    restart == 0 && (restart = n)

    GMRESSolver(op, maxiter, restart, tol)
end


operator(solver::GMRESSolver) = solver.linear_operator


function solve(solver::GMRESSolver, b)
    op = operator(solver)
    x, ch = IterativeSolvers.gmres(op, b, log=true,  maxiter=solver.maxiter,
        restart=solver.restart, tol=solver.tol)
    return x, ch
end


function Base.:*(solver::GMRESSolver, b)
    x, ch = solve(solver, b)
    println("Number of iterations: ", ch.iters)
    ch.isconverged || error("Iterative solver did not converge.")
    return x
end


function gmres(eq::DiscreteEquation; maxiter=0, restart=0, tol=0)

    test_space_dict  = eq.test_space_dict
    trial_space_dict = eq.trial_space_dict

    lhs = eq.equation.lhs
    rhs = eq.equation.rhs

    b = assemble(rhs, test_space_dict)
    Z = assemble(lhs, test_space_dict, trial_space_dict)

    if tol == 0
        invZ = GMRESSolver(Z, maxiter=maxiter, restart=restart)
    else
        invZ = GMRESSolver(Z, maxiter=maxiter, restart=restart, tol=tol)
    end
    x = invZ * b
end
