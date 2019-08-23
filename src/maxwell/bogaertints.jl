function momintegrals!(op::MWSingleLayer3D, g::RTRefSpace, f::RTRefSpace, t, s, z, strat::BogaertStrategy)

    T, GG  = GetIntegrals(t, s, op.gamma, strat)

    # Get the (u,v,w) → (x,y,z) tf matrix for tcell == bcell
    P = @SMatrix [
        t[1][1] t[2][1] t[3][1]
        t[1][2] t[2][2] t[3][2]
        t[1][3] t[2][3] t[3][3]
    ]

    R = @SMatrix [
        s[1][1] s[2][1] s[3][1]
        s[1][2] s[2][2] s[3][2]
        s[1][3] s[2][3] s[3][3]
    ]

    ∫G = sum(T)

    Q = P*T
    ∫xG = @SVector [
        Q[1,1]+Q[1,2]+Q[1,3],
        Q[2,1]+Q[2,2]+Q[2,3],
        Q[3,1]+Q[3,2]+Q[3,3],
    ]

    Q = R*transpose(T)
    ∫Gy = @SVector [
        Q[1,1]+Q[1,2]+Q[1,3],
        Q[2,1]+Q[2,2]+Q[2,3],
        Q[3,1]+Q[3,2]+Q[3,3]
    ]

    Q = P*T*transpose(R)
    ∫xGy = Q[1,1] + Q[2,2] + Q[3,3]

    c₁ = op.α
    c₂ = op.β

    # Now build the shape-shape interactions from these
    α = 1 / volume(t) / volume(s) / 4
    γ = op.gamma
    for i in 1:3
        a = t[i]
        for j in 1:3
            b = s[j]
            z[i,j] = α*c₁*(∫xGy - dot(a,∫Gy) - dot(b,∫xG) + dot(a,b)*∫G) + 4α*c₂*∫G
        end
    end

end



function momintegrals!(op::MWDoubleLayer3D, g::RTRefSpace, f::RTRefSpace,
    τ, σ, z, strat::BogaertStrategy)

    # Get the primitives
    r = τ.vertices
    s = σ.vertices

    G, GG = GetIntegrals(τ, σ, op.gamma, strat)

    # representation of RT elements in terms of
    # the tangents tu, tv and coordinates u,v,w
    # f_1 = [ (-v-w) t_u + v t_v ] / j
    # f_2 = [ u t_u + (-u-w) t_v ] / j
    # f_3 = [ u t_u + v t_v ] / j

    α = @SMatrix [
        0 -1 -1
        1  0  0
        1  0  0]

    β = @SMatrix [
         0  1  0
        -1  0 -1
         0  1  0]

    duu = τ.tangents[1] × σ.tangents[1]
    duv = τ.tangents[1] × σ.tangents[2]
    dvu = τ.tangents[2] × σ.tangents[1]
    dvv = τ.tangents[2] × σ.tangents[2]

    J = 4 * volume(σ) * volume(τ)
    for i in 1:3
        for j in 1:3
            Iuu = sum([α[i,p] * GG[p,q] * α[j,q] for p in 1:3, q in 1:3])
            Iuv = sum([α[i,p] * GG[p,q] * β[j,q] for p in 1:3, q in 1:3])
            Ivu = sum([β[i,p] * GG[p,q] * α[j,q] for p in 1:3, q in 1:3])
            Ivv = sum([β[i,p] * GG[p,q] * β[j,q] for p in 1:3, q in 1:3])
            z[i,j] = (dot(duu, Iuu) + dot(duv, Iuv) + dot(dvu, Ivu) + dot(dvv, Ivv)) / J
        end
    end

end
