# T: coeff type
# Degree: degree
# Dim1: dimension of the support + 1
struct LagrangeRefSpace{T,Degree,Dim1,NF} <: RefSpace{T} end

numfunctions(s::LagrangeRefSpace{T,D,2}, ch::CompScienceMeshes.ReferenceSimplex{1}) where {T,D} = D+1
numfunctions(s::LagrangeRefSpace{T,0,3}, ch::CompScienceMeshes.ReferenceSimplex{2}) where {T} = 1
numfunctions(s::LagrangeRefSpace{T,1,3}, ch::CompScienceMeshes.ReferenceSimplex{2}) where {T} = 3
numfunctions(s::LagrangeRefSpace{T,2,3}, ch::CompScienceMeshes.ReferenceSimplex{2}) where {T} = 6
numfunctions(s::LagrangeRefSpace{T,Dg}, ch::CompScienceMeshes.ReferenceSimplex{D}) where {T,Dg,D} = binomial(D+Dg,Dg)

# valuetype(ref::LagrangeRefSpace{T}, charttype) where {T} =
#         SVector{numfunctions(ref), Tuple{T,T}}
valuetype(ref::LagrangeRefSpace{T}, charttype) where {T} = T

# Evaluate constant lagrange elements on anything
(ϕ::LagrangeRefSpace{T,0})(tp) where {T} = SVector(((value=one(T), derivative=zero(T)),))
(ϕ::LagrangeRefSpace{T,0,3})(tp) where {T} = SVector(((value=one(T), derivative=zero(T)),))


# Evaluate linear Lagrange elements on a segment
function (f::LagrangeRefSpace{T,1,2})(mp) where T
    u = mp.bary[1]
    j = jacobian(mp)
    SVector(
        (value=  u, derivative=-1/j),
        (value=1-u, derivative= 1/j))
end

# Evaluete linear lagrange elements on a triangle
function (f::LagrangeRefSpace{T,1,3})(t) where T
    u,v,w, = barycentric(t)

    j = jacobian(t)
    p = t.patch
    σ = sign(dot(normal(t), cross(p[1]-p[3],p[2]-p[3])))
    SVector(
        (value=u, curl=σ*(p[3]-p[2])/j),
        (value=v, curl=σ*(p[1]-p[3])/j),
        (value=w, curl=σ*(p[2]-p[1])/j))
end



# Evaluate constant Lagrange elements on a triangle, with their curls
function (f::LagrangeRefSpace{T,0,3})(t, ::Type{Val{:withcurl}}) where T
    i = one(T)
    z = zero(cartesian(t))
    SVector(((value=i, curl=z,),))
end


function curl(ref::LagrangeRefSpace{T,1,3} where {T}, sh, el)
    sh1 = Shape(sh.cellid, mod1(sh.refid+1,3), -sh.coeff)
    sh2 = Shape(sh.cellid, mod1(sh.refid+2,3), +sh.coeff)
    return [sh1, sh2]
end

function gradient(ref::LagrangeRefSpace{T,1,4}, sh, tet) where {T}
    this_vert = tet.vertices[sh.refid]
    # other_verts = deleteat(tet.vertices, sh.refid)
    # opp_face = simplex(other_verts...)
    opp_face = faces(tet)[sh.refid]
    ctr_opp_face = center(opp_face)
    n = normal(ctr_opp_face)
    h = -dot(this_vert - cartesian(ctr_opp_face), n)
    @assert h > 0
    gradval = -(1/h)*n
    output = Vector{Shape{T}}()
    for (i,edge) in enumerate(CompScienceMeshes.edges(tet))
        ctr_edge = center(edge)
        tgt = tangents(ctr_edge,1)
        tgt = normalize(tgt)
        lgt = volume(edge)
        cff = -lgt * dot(tgt, gradval)
        isapprox(cff, 0, atol=sqrt(eps(T))) && continue
        push!(output, Shape(sh.cellid, i, sh.coeff * cff))
    end
    return output
end

function gradient(ref::LagrangeRefSpace{T,1,3} where {T}, sh, el)
    sh1 = Shape(sh.cellid, mod1(sh.refid+1,3), +sh.coeff)
    sh2 = Shape(sh.cellid, mod1(sh.refid+2,3), -sh.coeff)
    return [sh1, sh2]
end

function gradient(ref::LagrangeRefSpace{T,1,2}, sh, seg) where {T}

    sh.refid == 1 && return [Shape(sh.cellid, 1, +sh.coeff/volume(seg))]
    @assert sh.refid == 2
    return [Shape(sh.cellid, 1, -sh.coeff/volume(seg))]

end



function strace(x::LagrangeRefSpace, cell, localid, face)

    Q = zeros(scalartype(x),2,3)

    p1 = neighborhood(face, 1)
    p2 = neighborhood(face, 0)

    u1 = carttobary(cell, cartesian(p1))
    u2 = carttobary(cell, cartesian(p2))

    P1 = neighborhood(cell, u1)
    P2 = neighborhood(cell, u2)

    vals1 = x(P1)
    vals2 = x(P2)

    num_shapes = numfunctions(x, domain(cell))
    for j in 1:num_shapes
        Q[1,j] = vals1[j].value
        Q[2,j] = vals2[j].value
    end

    Q
end

function strace(x::LagrangeRefSpace{T, 1, 4, 4}, cell, localid, face) where {T}

    #T = scalartype(x)
    t = zeros(T, 3, 4)
    for (k,fvert) in enumerate(face.vertices)
        for (l,cvert) in enumerate(cell.vertices)
            nrm = norm(fvert - cvert)
            if isapprox(nrm, 0, atol=sqrt(eps(T)))
                t[k,l] = T(1.0)
                break
            end
        end
    end

    return t
end


function restrict(refs::LagrangeRefSpace{T,0}, dom1, dom2) where T
    n = numfunctions(refs, domain(dom1))
    Q = Matrix{T}(I, n, n)
end

function restrict(f::LagrangeRefSpace{T,1}, dom1, dom2) where T

    D = numfunctions(f, domain(dom1))
    Q = zeros(T, D, D)

    # for each point of the new domain
    for i in 1:D
        v = dom2.vertices[i]

        # find the barycentric coordinates in dom1
        uvn = carttobary(dom1, v)

        # evaluate the shape functions in this point
        x = neighborhood(dom1, uvn)
        fx = f(x)

        for j in 1:D
            Q[j,i] = fx[j][1]
        end
    end

    return Q
end



## Quadratic Lagrange element on a triangle
function (f::LagrangeRefSpace{T,2,3})(t) where T
    u,v,w, = barycentric(t)

    j = jacobian(t)
    p = t.patch

    σ = sign(dot(normal(t), cross(p[1]-p[3],p[2]-p[3])))
     SVector(
        (value=u*(2*u-1), curl=σ*(p[3]-p[2])*(4u-1)/j),
        (value=v*(2*v-1), curl=σ*(p[1]-p[3])*(4v-1)/j),
        (value=w*(2*w-1), curl=σ*(p[2]-p[1])*(4w-1)/j),
        (value=4*v*w, curl=4*σ*(w*(p[1]-p[3])+v*(p[2]-p[1]))/j),
        (value=4*w*u, curl=4*σ*(w*(p[3]-p[2])+u*(p[2]-p[1]))/j),
        (value=4*u*v, curl=4*σ*(u*(p[1]-p[3])+v*(p[3]-p[2]))/j),
    )
end


function curl(ref::LagrangeRefSpace{T,2,3} where {T}, sh, el)
    #curl of lagc0d2 as combination of bdm functions 
    z=zero(typeof(sh.coeff))
    if sh.refid < 4
        sh1 = Shape(sh.cellid, mod1(2*sh.refid+1,6), +sh.coeff)
        sh2 = Shape(sh.cellid, mod1(2*sh.refid+2,6), -3*sh.coeff)
        sh3 = Shape(sh.cellid, mod1(2*sh.refid+3,6), +3*sh.coeff)
        sh4 = Shape(sh.cellid, mod1(2*sh.refid+4,6), -sh.coeff)
    else
        sh1 = Shape(sh.cellid, mod1(2*sh.refid+4,6), z*sh.coeff)
        sh2 = Shape(sh.cellid, mod1(2*sh.refid+5,6), -4*sh.coeff)
        sh3 = Shape(sh.cellid, mod1(2*sh.refid+6,6), +4*sh.coeff)
        sh4 = Shape(sh.cellid, mod1(2*sh.refid+7,6), z*sh.coeff)
    end
    return [sh1, sh2, sh3, sh4]
end


const _vert_perms_lag = [
    (1,2,3),
    (2,3,1),
    (3,1,2),
    (2,1,3),
    (1,3,2),
    (3,2,1),
]

const _dof_perms_lag0 = [
    (1),
    (1),
    (1),
    (1),
    (1),
    (1),
]
const _dof_perms_lag1 = [
    (1,2,3),
    (3,1,2),
    (2,3,1),
    (2,1,3),
    (1,3,2),
    (3,2,1),
]

function dof_permutation(::LagrangeRefSpace{<:Any,0}, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_lag)
    return _dof_perms_lag0[i]
end

function dof_permutation(::LagrangeRefSpace{<:Any,1}, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_lag)
    return _dof_perms_lag1[i]
end

function dof_perm_matrix(::LagrangeRefSpace{<:Any,0}, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_rt)
    @assert i != nothing
    return _dof_lag0perm_matrix[i]
end

function dof_perm_matrix(::LagrangeRefSpace{<:Any,1}, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_rt)
    @assert i != nothing
    return _dof_rtperm_matrix[i]
end

const _dof_lag0perm_matrix = [
    @SMatrix[1],         # 1. {1,2,3}
    @SMatrix[1],         # 2. {2,3,1}
    @SMatrix[1],         # 3. {3,1,2}
    @SMatrix[1],         # 4. {2,1,3}
    @SMatrix[1],         # 5. {1,3,2}
    @SMatrix[1]         # 6. {3,2,1}
]

function (ϕ::LagrangeRefSpace{T, 1, 4, 4})(lag) where T

    u, v, w = parametric(lag)

    tu = tangents(lag, 1)
    tv = tangents(lag, 2)
    tw = tangents(lag, 3)

    B = [tu tv tw]
    A = inv(transpose(B))

    # gradient in u,v,w (unit tetrahedron)
    gr1=SVector{3, T}(1.0, 0.0, 0.0)
    gr2=SVector{3, T}(0.0, 1.0, 0.0)
    gr3=SVector{3, T}(0.0, 0.0, 1.0)
    gr4=SVector{3, T}(-1.0, -1.0, -1.0)

    return SVector((
        (value = u, gradient = A*gr1),
        (value = v, gradient = A*gr2),
        (value = w, gradient = A*gr3),
        (value = T(1.0)-u-v-w, gradient = A*gr4)
    ))
end



# Evaluate higher order Lagrange elements on triangles
# TODO: Optimise using code generation
function (ϕ::LagrangeRefSpace{T,Degree,3})(p) where {T,Degree}

    u, v = parametric(p)
    w = 1 - u - v
    idx = 0

    suppdim = 2
    localdim = binomial(suppdim+Degree, suppdim)
    vals = T[]
    diffus = T[]
    diffvs = T[]

    D1 = Degree + 1
    s = range(zero(T), one(T), length=D1)
    for i in 0:Degree
        ui = i/Degree
        for j in 0:Degree
            vj = j/Degree
            for k in 0:Degree
                wk = k/Degree
                i + j + k == Degree || continue

                prod_p = one(T)
                for p in 0:i-1
                    up = p / Degree
                    prod_p *= (u-up) / (ui-up)
                end
                prod_q = one(T)
                for q in 0:j-1
                    vq = q / Degree
                    prod_q *= (v-vq) / (vj-vq)
                end
                prod_r = one(T)
                for r in 0:k-1
                    wr = r / Degree
                    prod_r *= (w-wr) / (wk-wr)
                end
                push!(vals, prod_p * prod_q * prod_r)

                diffu = zero(T)
                diffv = zero(T)
                for l in 0:i-1
                    ul = l/Degree
                    prod_pl = one(T)
                    for p in 0:i-1
                        p == l && continue
                        up = p/Degree
                        prod_pl *= (u-up) / (ui-up)
                    end
                    diffu += prod_pl * prod_q * prod_r / (ui-ul)
                end
                for m in 0:j-1
                    vm = m/Degree
                    prod_qm = one(T)
                    for q in 0:j-1
                        q == m && continue
                        vq = q/Degree
                        prod_qm *= (v-vq) / (vj-vq)
                    end
                    diffv += prod_p * prod_qm * prod_r / (vj-vm)
                end
                for n in 0:k-1
                    wn = n/Degree
                    prod_rn = one(T)
                    for r in 0:k-1
                        r == n && continue
                        wr = r/Degree
                        prod_rn *= (w-wr) / (wk-wr)
                    end
                    diffu -= prod_p * prod_q * prod_rn / (wk-wn)
                    diffv -= prod_p * prod_q * prod_rn / (wk-wn)
                end

                push!(diffus, diffu)
                push!(diffvs, diffv)

                idx += 1
    end end end
 
    tu = tangents(p,1)
    tv = tangents(p,2)
    j = jacobian(p)
    NF = length(vals)
    SVector{NF}([(value=f, curl=(-dv*tu+du*tv)/j) for (f,du,dv) in zip(vals, diffus, diffvs)])
end

# fields[i] ≈ sum(Q[j,i] * interpolant[j].value for j in 1:numfunctions(interpolant))
function interpolate(fields, interpolant::LagrangeRefSpace{T,Degree,3}, chart) where {T,Degree}

    dim = binomial(2+Degree, Degree)

    I = 0:Degree
    s = range(0,1,length=Degree+1)
    Is = zip(I,s)
    idx = 1
    vals = Vector{Vector{T}}()
    for (i,ui) in Is
        for (j,vj) in Is
            for (k,wk) in Is
                i + j + k == Degree || continue
                @assert ui + vj + wk ≈ 1
                p = neighborhood(chart, (ui,vj))
                push!(vals, fields(p))
                idx += 1
    end end end

    # Q = hcat(vals...)
    Q = Matrix{T}(undef, length(vals[1]), length(vals))
    for i in eachindex(vals)
        Q[:,i] .= vals[i]
    end
    return Q
end