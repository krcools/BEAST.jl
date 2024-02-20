function update!(op::ConductivityTD_Functionaltype, jcoeffs, ecoeffs, k, eq1, eq2)
    tgeo = BEAST.geometry(eq2.test_space_dict[1].space)
    bgeo = BEAST.geometry(eq2.trial_space_dict[1].space)
    if CompScienceMeshes.refines(tgeo, bgeo)
        op.efield = quadpoint_field_refines(op, ecoeffs, k, eq2.test_space_dict[1], eq2.trial_space_dict[1])
        op.jflux = quadpoint_field_refines(op, jcoeffs, k, eq2.test_space_dict[1], eq1.trial_space_dict[1])
    else
        op.efield = quadpoint_field(op, ecoeffs, k, eq2.test_space_dict[1], eq2.trial_space_dict[1])
        op.jflux = quadpoint_field(op, jcoeffs, k, eq2.test_space_dict[1], eq1.trial_space_dict[1])
    end	
end

function mpeval(point, op::ConductivityTD_Functionaltype; type=nothing)
    tcell = op.tcells[op.timeindex]
    k = op.timeindex
    #i = CompScienceMeshes.findchart(op.charts, op.chart_tree, point)
    i=nothing
    if k==1
        i = CompScienceMeshes.findchart(op.charts, op.chart_tree, point)
        push!(op.ptlocs, i)
        push!(op.pts, cartesian(point))
    else
        op.lastloc +=1
        if isapprox(op.pts[op.lastloc],cartesian(point), atol=1e-12)
            i=op.ptlocs[op.lastloc]
        else
            println("incorrect point related from last iteration")
            i = CompScienceMeshes.findchart(op.charts, op.chart_tree, point)
        end
    end
    if i !== nothing
        # @show i
        chart = op.charts[i]
        u = carttobary(chart, point)
        @assert isapprox(cartesian(point), barytocart(chart, u), atol=1e-12)
        vals = op.refs(neighborhood(chart,u))
        tvals = op.trefs(neighborhood(tcell,1))
        op.field[1] = [0.0,0.0,0.0]
        for si in 1:op.numrefs
            for ti in 1:op.tnumrefs
                for (m,a) in op.ad[i,si]
                    for (n,b) in op.tad[k,ti]
                        fx = vals[si].value
                        tfx = tvals[ti]
                        op.field[1] += (op.coeffs[m,n] * tfx * b) * a * fx
                    end
                end
            end
        end
    end
    return op.field[1]
end

function quadpoint_field_refines(biop::ConductivityTD_Functionaltype, coeffs, k, tfs, bfs;
        quadstrat=BEAST.defaultquadstrat(biop, tfs.space, bfs.space))
    
        tgeo = BEAST.geometry(tfs.space)
        bgeo = BEAST.geometry(bfs.space)
        @assert CompScienceMeshes.refines(tgeo, bgeo)
    
        trefs = BEAST.refspace(tfs.space)
        brefs = BEAST.refspace(bfs.space)
        btrefs = BEAST.refspace(bfs.time)

        numbrefs = BEAST.numfunctions(brefs)
        numbtrefs = BEAST.numfunctions(btrefs)
    
        tels, tad, ta2g = BEAST.assemblydata(tfs.space)
        bels, bad, ba2g = BEAST.assemblydata(bfs.space)
        btels, btad = BEAST.assemblydata(bfs.time)

        bg2a = zeros(Int, length(BEAST.geometry(bfs.space)))
        for (i,j) in enumerate(ba2g) bg2a[j] = i end
    
        qd = BEAST.quaddata(biop, trefs, brefs, tels, bels, quadstrat)
        T = eltype(coeffs)
        D = dimension(BEAST.geometry(bfs.space))
        U = D+1
        PT = SVector{U, T}
        field = zeros(PT, (length(tels),length(qd[1])))
        btcell = btels[k]
        for (p,tcell) in enumerate(tels)
    
            P = ta2g[p]
            Q = CompScienceMeshes.parent(tgeo, P)
            q = bg2a[Q]
    
            bcell = bels[q]
            @assert overlap(tcell, bcell)
    
            isct = intersection(tcell, bcell)
            for cell in isct    
                qr = BEAST.quadrule(biop, trefs, brefs, cell, qd, quadstrat)
                for (qi,qdpt) in enumerate(qr)
                    mp = carttobary(bcell, qdpt[2])
                    vals = brefs(neighborhood(bcell, mp))
                    tvals = btrefs(neighborhood(btcell,1))
                    for si in 1:numbrefs
                        fx = vals[si].value
                        for ti in 1:numbtrefs
                            tfx = tvals[ti]
                            for (m,a) in bad[q,si]
                                for (n,b) in btad[k,ti]
                                    field[p,qi] += (coeffs[m,n] * tfx * b) * a * fx
                                end
                            end
                        end
                    end
                end
            end # next cell in intersection
        end # next cell in the parent geometry
        return field
end

function quadpoint_field(biop::ConductivityTD_Functionaltype, coeffs, k, tfs, bfs;
    quadstrat=BEAST.defaultquadstrat(biop, tfs.space, bfs.space))

    trefs = BEAST.refspace(tfs.space)
    brefs = BEAST.refspace(bfs.space)
    btrefs = BEAST.refspace(bfs.time)

    numbrefs = BEAST.numfunctions(brefs)
    numbtrefs = BEAST.numfunctions(btrefs)

    tels, tad, ta2g = BEAST.assemblydata(tfs.space)
    bels, bad, ba2g = BEAST.assemblydata(bfs.space)
    btels, btad = BEAST.assemblydata(bfs.time)

    bg2a = zeros(Int, length(BEAST.geometry(bfs.space)))
    for (i,j) in enumerate(ba2g) bg2a[j] = i end

    qd = BEAST.quaddata(biop, trefs, brefs, tels, bels, quadstrat)
    T = eltype(coeffs)
    D = dimension(BEAST.geometry(bfs.space))
    U = D+1
    PT = SVector{U, T}
    field = zeros(PT, (length(tels),length(qd[1])))
    btcell = btels[k]
    for (p,bcell) in enumerate(bels)   
        qr = BEAST.quadrule(biop, trefs, brefs, bcell, qd, quadstrat)
        #println(length(qr))
        for (qi,qdpt) in enumerate(qr)
            mp = carttobary(bcell, qdpt[2])
            vals = brefs(neighborhood(bcell, mp))
            tvals = btrefs(neighborhood(btcell,1))
            for si in 1:numbrefs
                fx = vals[si].value
                for ti in 1:numbtrefs
                    tfx = tvals[ti]
                    for (m,a) in bad[p,si]
                        for (n,b) in btad[k,ti]
                            field[p,qi] += (coeffs[m,n] * tfx * b) * a * fx
                        end
                    end
                end
            end
        end
    end
    return field
end

function marchonintimenl(eq1, eq2,  Z, inc, Ġ, G_j, G_nl, Nt)
    T = eltype(Z)
    Z0 = zeros(T, size(Z)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Z0,Z,1)
    Ġ0 = Ġ.data[1,:,:]
    G_j0 = G_j.data[1,:,:]
    G_nl0 = G_nl.data[1,:,:]
    T = eltype(Z0)
    M,N = size(Z0)
    Me,Ne = size(Ġ0)
    V0 = zeros(N+Ne,N+Ne)
    V0[1:N, 1:N] = Z0
    V0[1:N, N+1:N+Ne] = -Ġ0
    V0[N+1:N+Ne, 1:N] = G_j0
    sol = zeros(2*N)
    xj_all = zeros(N)
    xe_all = zeros(N)
    @assert M == size(inc,1)
    xj = zeros(T,N,Nt)
	xe = zeros(T,Ne,Nt)
	xeprev = zeros(T,N)
    yj = zeros(T,N)
	ye = zeros(T,N)
    jk_start = 2
    ek_start = 2
    jk_stop = Nt
    ek_stop = BEAST.numfunctions(eq1.trial_space_dict[1].time)+1
    csxj = zeros(T,N,Nt)
    csxe = zeros(T,Ne,Nt)
    σ = eq2.equation.rhs.terms[1].functional
    σop = BEAST.ConductivityTDOp(σ)
    bσ = zeros(T, N)
    iZ = BEAST.GMRESSolver(Z0, restart=0, reltol=1e-6)
    invZ = inv(Z0)
    C1 = G_j0*invZ
    C2 = C1*Ġ0
    #= try =#
    for i in 1:Nt
        println(i)
        R = inc[:,i]
        fill!(yj,0)
        fill!(ye,0)
        BEAST.ConvolutionOperators.convolve!(yj,Z,xj,csxj,i,jk_start,jk_stop)
        BEAST.ConvolutionOperators.convolve!(ye,Ġ,xe,csxe,i,ek_start,ek_stop)
        iter_max = 100
        for l in 1:iter_max
            if l==1
                if i==1
                    i=2
                end
                xeprev = xe[:,i-1]
                update!(σ, xj, xe, i-1, eq1, eq2)
            else
                xeprev = xe[:,i]
                update!(σ, xj, xe, i, eq1, eq2)
            end
            bσ = BEAST.assemble(σ, eq2.test_space_dict[1].space)
            Q = BEAST.assemble(σop, eq2.test_space_dict[1].space, eq2.trial_space_dict[1].space)
            rhs1 = R - yj + ye
            rhs2 = bσ - Q*xeprev
            V0[N+1:N+Ne,N+1:N+Ne] = -Q
            b = [rhs1; rhs2]
            sol = inv(Matrix(V0))*b
            xj[:,i] = sol[1:N]
            xe[:,i] = sol[N+1:end]
            #= rhs1 = R - yj + ye
            rhs2 = bσ - Q*xeprev - C1*rhs1
            iw = BEAST.GMRESSolver(C2-Q, restart=0, reltol=1e-6)
            u, ch = solve(iw, rhs2)
            xe[:,i] .= u
            xj[:,i] .= invZ*(rhs1+Ġ0*xe[:,i]) =#
            #xe_all = hcat(xe_all, xe[:,i])
            #xj_all = hcat(xj_all, xj[:,i])
            println("norm xe ", norm(xe[:,i]-xeprev)/M)
            if norm(xe[:,i]-xeprev)/M < 1e-4
                break
            end
            #= println("L_inf norm of diff of xe ", maximum(abs.((xe[:,i]-xeprev))))
            if maximum(abs.((xe[:,i]-xeprev)))/M < 1e-10
                break
            end =#
        end
        if i > 1
            csxj[:,i] .= csxj[:,i-1] .+ xj[:,i]
            csxe[:,i] .= csxe[:,i-1] .+ xe[:,i]
        else
            csxj[:,i] .= xj[:,i]
            csxe[:,i] .= xe[:,i]
        end
        (i % 10 == 0) && print(i, "[", Nt, "] - ")
    end
    return xj, xe, xj_all, xe_all
    #= catch e
        return xj, xe, xj_all, xe_all
    end =#
end

function td_solve(eq1::BEAST.DiscreteEquation, eq2::BEAST.DiscreteEquation)
    Z = BEAST.td_assemble(eq1.equation.lhs, eq1.test_space_dict, eq1.trial_space_dict)
    b = BEAST.td_assemble(eq1.equation.rhs, eq1.test_space_dict)
    h = eq2.trial_space_dict[eq2.equation.lhs.terms[1].trial_id]
    h1 = h.space⊗BEAST.derive(h.time)
    f1 = eq1.test_space_dict[1]
    f2 = eq2.test_space_dict[1]
    g = eq1.trial_space_dict[1]
    idST = BEAST.Identity()⊗BEAST.Identity()
    Ġ = BEAST.assemble(idST, f1, h1)
    G_j = BEAST.assemble(idST, f2, g)
    Nt = BEAST.numfunctions(g.time)
    if typeof(eq2.equation.rhs.terms[1].functional)==ConductivityTDFunc
        G_nl = BEAST.assemble(idST, f2, h)
        return marchonintimenl(eq1, eq2, Z, b, Ġ, G_j, G_nl, Nt)
    end
end