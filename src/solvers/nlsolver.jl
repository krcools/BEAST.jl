function newtoniterator(eq1, eq2, inc, Z, G_e, G_j)
    T = eltype(G_e)
    M,N = size(G_e)
    V0 = zeros(N+N,N+N)
    V0[1:N, 1:N] = 2.0*Z
    V0[1:N, N+1:N+N] = G_e 
    V0[N+1:N+N, 1:N] = G_j
    sol = zeros(2*N)
    @assert M == size(inc,1)
    xj = zeros(T,N,5)
	xe = zeros(T,N,5)
	xeprev = zeros(T,N,1)
    xj_all = zeros(N)
    xe_all = zeros(N)
    σ = eq2.equation.rhs.terms[1].functional
    σop = BEAST.ConductivityTDOp(σ)
    bσ = zeros(T, N)
    iter_max = 100
    R = inc[:,Int(floor(end/2))]
    for l in 1:iter_max
        xeprev = xe[:,1]
        update!(σ, xj, xe, 1, eq1, eq2)
        bσ = BEAST.assemble(σ, eq2.test_space_dict[1].space)
        Q = BEAST.assemble(σop, eq2.test_space_dict[1].space, eq2.trial_space_dict[1].space)
        #Q_inv = inv(Matrix(Q))
        #= println("normQ ", norm(Q))
        println("normbsigma ", norm(bσ)) =#
        #V0 = inv(Z0*G_j0_inv*Q - Ġ0)
        V0[N+1:N+N,N+1:N+N] = -Q
        rhs1 = R
        rhs2 = bσ - Q*xeprev
        #b = rhs1 + rhs2         #solve |Z   -Ġ|J =|b|
        b = [rhs1; rhs2]
        sol = inv(Matrix(V0))*b
        #= invV0 = GMRESSolver(V0, maxiter=1000, restart=0, tol=1e-6)
        mul!(sol, invV0, b) =#
        xj[:,1] = sol[1:N]                            #|G_j  Q|E =|0|
        xe[:,1] = sol[N+1:N+N]
        xe_all = hcat(xe_all, xe[:,1])
        xj_all = hcat(xj_all, xj[:,1])
        #= println(l, " norm xe ", norm(xe[:,1]-xeprev))
        if norm(xe[:,1]-xeprev) < 1e-6
            break
        end =#
        println("L_inf norm of diff of xe ", maximum(abs.((xe[:,1]-xeprev))))
        if maximum(abs.((xe[:,1]-xeprev))) < 1e-12
            #test if the discrete ohms law is satisfied
            update!(σ, xj, xe, 1, eq1, eq2)
            bσ = BEAST.assemble(σ, eq2.test_space_dict[1].space)
            println(norm(G_j*xj[:,1] - bσ))
            break
        end
    end
    return xj, xe, xj_all, xe_all
end

function nlsolve(eq1::BEAST.DiscreteEquation, eq2::BEAST.DiscreteEquation)
    #Z = BEAST.td_assemble(eq1.equation.lhs, eq1.test_space_dict, eq1.trial_space_dict)
    b = BEAST.td_assemble(eq1.equation.rhs, eq1.test_space_dict)
    h = eq2.trial_space_dict[eq2.equation.lhs.terms[1].trial_id]
    h1 = h.space
    f1 = eq1.test_space_dict[1].space
    f2 = eq2.test_space_dict[1].space
    g = eq1.trial_space_dict[1].space
    id = BEAST.Identity()
    Z = BEAST.assemble(id, f1, g)
    G_e = BEAST.assemble(id, f1, h1)
    G_j = BEAST.assemble(id, f2, g)
    return newtoniterator(eq1, eq2, b, Z, G_e, G_j)
end