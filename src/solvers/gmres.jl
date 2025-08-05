import Krylov

struct GMRES{T} <: LinearMap{T}
    A
    M
    N
    memory
    restart
    atol
    rtol
    itmax
    verbose
    workspace
end

Base.axes(solver::GMRES) = reverse(axes(solver.A))
Base.size(solver::GMRES) = reverse(size(solver.A))


function GMRES(A;
    M = LinearAlgebra.I,
    N = LinearAlgebra.I,
    memory = 20,
    restart = false,
    atol = sqrt(eps(real(eltype(A)))),
    rtol = sqrt(eps(real(eltype(A)))),
    itmax = 2 * size(A,2),
    verbose = 1)

    m, n = size(A)
    T = eltype(A)
    ws = Krylov.GmresWorkspace(m, n, Vector{T}; memory)
    return GMRES{T}(A, M, N, memory, restart, atol, rtol, itmax, verbose, ws)
end


function solve!(x, solver::GMRES, b)
    Krylov.gmres!(solver.workspace, solver.A, Vector(b);
        M = solver.M, 
        N = solver.N,
        # memory = solver.memory, 
        restart = solver.restart, 
        atol = solver.atol, 
        rtol = solver.rtol, 
        itmax = solver.itmax, 
        verbose = solver.verbose,
        history = true)
    x .= Krylov.solution(solver.workspace)
    return x, Krylov.statistics(solver.workspace)
end

function solve(solver::GMRES, b)
    T = promote_type(eltype(solver), eltype(b))
    x = similar(Vector{T}, size(solver)[2])
    fill!(x,0)
    x, ch = solve!(x, solver, b)
    # z = similar(Array{T}, axes(solver)[2])
    # copyto!(z,x)
    return x, ch
end

@testitem "GMRES" begin
    using LinearAlgebra
    
    A = rand(10, 10)
    b = rand(10)

    Ai = BEAST.GMRES(A; restart=true, atol=1e-8, rtol=1e-8, verbose=1, memory=200)
    x, ch = solve(Ai, b)

    Ai2 = BEAST.GMRESSolver(A; restart=200, abstol=1e-8, reltol=1e-8, verbose=true)
    x2, ch2 = solve(Ai2, b)

    @test norm(x - x2) < 1e-6
end