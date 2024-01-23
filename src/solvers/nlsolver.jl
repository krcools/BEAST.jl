function newtoniterator(eq1, eq2, inc, Z, G_e, G_j)
    T = eltype(G_j)
    M,N = size(G_e)
    #= Z0 = zeros(T, size(Z)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Z0,Z,1) =#
    #V0 = zeros(M+N,M+N)
    Z0 = zeros(T, size(Z))
    Z0 = 1.0*Z
    #V0[1:M, 1:M] = Z0
    #V0[1:M, M+1:M+N] = G_e 
    #V0[M+1:M+N, 1:M] = G_j
    sol = zeros(M+N)
    @assert M == size(inc,1)
    xj = zeros(T,M,5)
	xe = zeros(T,N,5)
	xeprev = zeros(T,N,1)
    xj_all = zeros(M)
    xe_all = zeros(N)
    σ = eq2.equation.rhs.terms[1].functional
    σop = BEAST.ConductivityTDOp(σ)
    bσ = zeros(T, N)
    iter_max = 100
    R = inc[:,50]
    invZ = inv(Matrix(Z0))
    C1 = G_j*invZ
    C2 = C1*G_e
    for l in 1:iter_max
        xeprev = xe[:,1]
        update!(σ, xj, xe, 1, eq1, eq2)
        bσ = BEAST.assemble(σ, eq2.test_space_dict[1].space)
        Q = BEAST.assemble(σop, eq2.test_space_dict[1].space, eq2.trial_space_dict[1].space)
        #Q_inv = inv(Matrix(Q))
        #= println("normQ ", norm(Q))
        println("normbsigma ", norm(bσ)) =#
        #V0 = inv(Z0*G_j0_inv*Q - Ġ0)
        V0 = [Z0 G_e; G_j -Q]
        #V0[M+1:M+N,M+1:M+N] = -Q
        rhs1 = R
        rhs2 = bσ - Q*xeprev
        #= rhs2 = bσ - Q*xeprev - C1*rhs1
        iw = BEAST.GMRESSolver(-C2-Q, restart=0, reltol=1e-6)
        u, ch = solve(iw, rhs2)
        xe[:,1] .= u
        xj[:,1] .= invZ*(rhs1-G_e*xe[:,1]) =#
        #b = rhs1 + rhs2         #solve |Z   -Ġ|J =|b|
        b = [rhs1; rhs2]
        #sol = inv(Matrix(V0))*b
        sol = V0\b
        #invV0 = BEAST.GMRESSolver(V0, maxiter=5000, restart=0, reltol=1e-6)
        #mul!(sol, invV0, b)
        xj[:,1] = sol[1:M]                            #|G_j  Q|E =|0|
        xe[:,1] = sol[M+1:M+N]
        xe_all = hcat(xe_all, xe[:,1])
        xj_all = hcat(xj_all, xj[:,1])
        #= println(l, " norm xe ", norm(xe[:,1]-xeprev))
        if norm(xe[:,1]-xeprev) < 1e-6
            break
        end =#
        println("L_inf norm of diff of xe ", maximum(abs.((xe[:,1]-xeprev))))
        if maximum(abs.((xe[:,1]-xeprev))) < 1e-4
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