import IterativeSolvers



struct GMRESSolver{T,L,R,PL,PR} <: LinearMap{T}
    linear_operator::L
    maxiter::Int
    restart::Int
    abstol::R
    reltol::R
    verbose::Bool
    left_preconditioner::PL
    right_preconditioner::PR
end


Base.axes(A::GMRESSolver) = reverse(axes(A.linear_operator))


function GMRESSolver(op::L;
    left_preconditioner = nothing,
    right_preconditioner = nothing,
    Pl = nothing,
    Pr = nothing,
    maxiter=0,
    restart=0,
    abstol::R = zero(real(eltype(op))),
    reltol::R = sqrt(eps(real(eltype(op)))),
    verbose=true) where {L,R<:Real}

    if left_preconditioner == nothing
        Pl == nothing && (Pl = IterativeSolvers.Identity())
    else
        if Pl == nothing
            Pl = BEAST.Preconditioner(left_preconditioner)
        else
            error("Either supply Pl or left_preconditioner, not both.")
        end
    end
    @assert Pl != nothing

    if right_preconditioner == nothing
        Pr == nothing && (Pr = IterativeSolvers.Identity())
    else
        if Pr == nothing
            Pr = BEAST.Preconditioner(right_preconditioner)
        else
            error("Either supply Pl or left_preconditioner, not both.")
        end
    end
    @assert Pr != nothing

    m, n = size(op)
    @assert m == n

    maxiter == 0 && (maxiter = n)
    restart == 0 && (restart = n)

    PL = typeof(Pl)
    PR = typeof(Pr)
    T = eltype(op)
    GMRESSolver{T,L,R,PL,PR}(op, maxiter, restart, abstol, reltol, verbose, Pl, Pr)
end


operator(solver::GMRESSolver) = solver.linear_operator


function solve(solver::GMRESSolver, b; abstol=solver.abstol, reltol=solver.reltol)
    T = promote_type(eltype(solver), eltype(b))
    x = similar(Array{T}, axes(solver)[2])
    fill!(x,0)
    x, ch = solve!(x, solver, b; abstol, reltol)
end


function solve!(x, solver::GMRESSolver, b; abstol=solver.abstol, reltol=solver.reltol)
    op = operator(solver)
    x, ch = IterativeSolvers.gmres!(x, op, b; 
        log=true, 
        maxiter=solver.maxiter,
        restart=solver.restart,
        reltol=reltol,
        abstol=abstol,
        verbose=solver.verbose,
        Pl=solver.left_preconditioner,
        Pr=solver.right_preconditioner)
    return x, ch
end


function Base.:*(A::GMRESSolver, b::AbstractVector)

    x, ch = solve(A, b)
    return x
    # T = promote_type(eltype(A), eltype(b))
    # y = BlockedVector{T}(undef, (axes(A,2),))
    # y = similar(Array{T}, axes(A,2))

    # mul!(y, A, b)
end

Base.size(solver::GMRESSolver) = reverse(size(solver.linear_operator))

function LinearAlgebra.mul!(y::AbstractVecOrMat, solver::GMRESSolver, x::AbstractVector)
    fill!(y,0)
    y, ch = solve!(y, solver, x)
    solver.verbose && println("Number of iterations: ", ch.iters)
    ch.isconverged || error("Iterative solver did not converge.")
    return y
end

LinearAlgebra.adjoint(A::GMRESSolver) = GMRESSolver(adjoint(A.linear_operator); maxiter=A.maxiter, restart=A.restart, abstol=A.abstol, reltol=A.reltol, verbose=A.verbose, Pl=A.right_preconditioner, Pr=A.left_preconditioner)
LinearAlgebra.transpose(A::GMRESSolver) = GMRESSolver(transpose(A.linear_operator); maxiter=A.maxiter, restart=A.restart, abstol=A.abstol, reltol=A.reltol, verbose=A.verbose, Pl=A.right_preconditioner, Pr=A.left_preconditioner)



function gmres_ch(eq::DiscreteEquation; maxiter=0, restart=0, tol=0, verbose=true)

    lhs = eq.equation.lhs
    rhs = eq.equation.rhs

    X = _spacedict_to_directproductspace(eq.test_space_dict)
    Y = _spacedict_to_directproductspace(eq.trial_space_dict)

    b = assemble(rhs, X)
    Z = assemble(lhs, X, Y)

    if tol == 0
        invZ = GMRESSolver(Z; maxiter, restart, verbose)
    else
        invZ = GMRESSolver(Z; maxiter, restart, reltol=tol, verbose)
    end
    x, ch = solve(invZ, b)
    # x = invZ * b

    ax = nestedrange(Y, 1, numfunctions)
    return BlockedVector(x, (ax,)), ch
end

gmres(eq::DiscreteEquation; maxiter=0, restart=0, tol=0) = gmres_ch(eq; maxiter, restart, tol)[1]

@testitem "GMRESSolver blocked" begin
    using LinearAlgebra
    using BEAST.BlockArrays

    n = 3
    J = Matrix{Float64}(I,n,n)
    Z = zeros(Float64,n,n)
    A = BlockedArray([Z J; -J Z], [n,n], [n,n])

    Ai = BEAST.GMRESSolver(A)
    b = collect(1:2n)

    x = Ai * b
    @test blocksize(x) == (2,)
    @test blocksizes(x) == [(n,), (n,)]

    @test x[1:n] ≈ -b[n+1:end]
    @test x[n+1:end] ≈ b[1:n]
end

struct CGSolver{L,M,T} <: LinearMap{T}
    A::L
    Pl::M
    abstol::T
    reltol::T
    maxiter::Int
    verbose::Bool
end


Base.size(solver::CGSolver) = reverse(size(solver.A))
Base.axes(solver::CGSolver) = reverse(axes(solver.A))
operator(solver::CGSolver) = solver.A


cg(A;
    Pl = IterativeSolvers.Identity(),
    abstol::Real = zero(real(eltype(A))),
    reltol::Real = sqrt(eps(real(eltype(A)))),
    maxiter::Int = size(A,2),
    verbose::Bool = false) = CGSolver(A, Pl, abstol, reltol, maxiter, verbose)


function solve(solver::CGSolver, b)
    T = promote_type(eltype(solver), eltype(b))
    x = similar(Vector{T}, size(solver)[2])
    fill!(x,0)
    x, ch = solve!(x, solver, b)
    z = similar(Array{T}, axes(solver)[2])
    copyto!(z,x)
    return z, ch
end


function solve!(x, solver::CGSolver, b)
    op = operator(solver)
    x, ch = IterativeSolvers.cg!(x, op, b;
        Pl = solver.Pl,
        abstol = solver.abstol,
        reltol = solver.reltol,
        maxiter = solver.maxiter,
        verbose = solver.verbose,
        log = true)
    return x, ch
end

function LinearAlgebra.mul!(y::AbstractVecOrMat, solver::CGSolver, x::AbstractVector)
    fill!(y,0)
    y, ch = solve!(y, solver, x)
    return y
end


struct Preconditioner{L,T} <: LinearMap{T}
    A::L
end

Preconditioner(A::L) where {L} = Preconditioner{L,eltype(A)}(A)

function LinearAlgebra.ldiv!(y, P::Preconditioner, x)
    mul!(y, P.A, x)
end

function LinearAlgebra.ldiv!(P::Preconditioner, b)
    c = deepcopy(b)
    mul!(b, P.A, c)
end

function Base.size(p::Preconditioner)
    reverse(size(p.A))
end


function solve(iA::AbstractMatrix, b)
    ch = nothing
    x = iA * b
    return x, ch
end 