struct GWPCurlSpace{T,M,P} <: Space{T}
    geo::M
    fns::Vector{Vector{Shape{T}}}
    pos::P
    degree::Int
end

function refspace(s::GWPCurlSpace{T}) where {T} GWPCurlRefSpace{T,s.degree}() end
function subset(s::S,I) where {S<:GWPCurlSpace} S(s.geo, s.fns[I], s.pos[I], s.degree) end

struct GWPGlobalEdgeDoFs order::Int end
struct GWPGlobalFaceDoFs order::Int end

# function globaldofs(edge_ch, face_ch, localspace, dof::GWPGlobalEdgeDoFs)
#     f = trace(edge_ch, face_ch, localspace)
#     T = coordtype(edge_ch)
#     ds = one(T) / (dof.order+2)
#     return stack(range(ds, step=ds, length=dof.order+1)) do s
#         p = neighborhood(edge_ch, s)
#         t = tangents(p, 1)
#         [dot(t,x.value) for x in f(s)]
#     end
# end

function globaldofs(global_ch, face_ch, localspace, dof::GWPGlobalEdgeDoFs)

    d = dof.order
    T = coordtype(global_ch)

    f = trace(global_ch, face_ch, localspace)
    g = u -> collect(fi.value for fi in f(u))
    Q = zeros(T, numfunctions(localspace, domain(face_ch)), d+1)
    return globaldofs!(Q, global_ch, g, dof)
end


function globaldofs!(Q, global_ch, pulled_back_fn, dof::GWPGlobalEdgeDoFs)

    d = dof.order
    T = coordtype(global_ch)
    ds = one(T) / (d+2)

    for i in 1:d+1
        s = i * ds
        p = neighborhood(global_ch, s)
        t = tangents(p, 1)
        vals = pulled_back_fn(s)
        for j in axes(Q,1)
            # Q[j,i] = dot(t, vals[j].value)
            # @show t, vals[j]
            Q[j,i] = dot(t, vals[j])
        end
    end

    return Q
end


@testitem "GWP edge global dofs" begin
    using CompScienceMeshes
    T = Float64
    face_ch = simplex(
        point(3,0,0),
        point(2,0,-1),
        point(0,0,1),
    )
    edge_ch = simplex(
        point(0,0,1),
        point(2,0,-1),
    )

    degree = 3
    localspace = BEAST.GWPCurlRefSpace{T,degree}()
    dofs = BEAST.GWPGlobalEdgeDoFs(degree)
    A = BEAST.globaldofs(edge_ch, face_ch, localspace, dofs)
    # display(round.(A, digits=3))
end

"""
    globaldofs(global_ch, face_ch, localspace, dof::GWPGlobalFaceDoFs)

Compute the face-interior global degrees of freedom for a local GWP curl space.
The returned matrix has one row per local basis function and two columns per
interior lattice point, containing the basis functions' components along the
two face tangents at that point.
"""
function globaldofs(global_ch, face_ch, localspace, dof::GWPGlobalFaceDoFs)

    d = dof.order
    T = coordtype(global_ch)

    f = trace(global_ch, face_ch, localspace)
    g = u -> collect(fi.value for fi in f(u))
    Q = zeros(T, numfunctions(localspace, domain(face_ch)), d*(d+1))
    return globaldofs!(Q, global_ch, g, dof)
end

function globaldofs!(Q, global_ch, pulled_back_fn, dof::GWPGlobalFaceDoFs)

    d = dof.order
    T = coordtype(global_ch)
    ds = one(T) / (d+2)

    S = ((i,j,d+2-i-j) for i in 1:d+1 for j in 1:d+1 if d+2-i-j > 0)
    idx = 1
    for (i,j,k) in S
        u = (i*ds,j*ds)
        p = neighborhood(global_ch, u)
        t1 = tangents(p,1)
        t2 = tangents(p,2)
        vals = pulled_back_fn(u)
        q1 = [dot(t1,x) for x in vals]
        q2 = [dot(t2,x) for x in vals]
        Q[:,idx] = q1; idx += 1
        Q[:,idx] = q2; idx += 1
    end

    return Q
end


@testitem "GWP face global dofs" begin
    using CompScienceMeshes
    using LinearAlgebra

    T = Float64
    face_ch = simplex(
        point(3,0,0),
        point(2,0,-1),
        point(0,0,1),
    )

    global_ch = simplex(
        point(3,0,0),
        point(2,0,-1),
        point(0,0,1),
    )

    degree = 3
    localspace = BEAST.GWPCurlRefSpace{T,degree}()
    dofs = BEAST.GWPGlobalFaceDoFs(degree)
    A = BEAST.globaldofs(global_ch, face_ch, localspace, dofs)
    num_inner_dofs = degree*(degree+1)
    T = round.(A[end-num_inner_dofs+1:end,:], digits=3)
    @test T ≈ Matrix{Int}(I, num_inner_dofs, num_inner_dofs)
end

function _addshapes!(fns, cell, gids, lids, β)
    T = eltype(β)
    atol = sqrt(eps(T))
    for i in axes(β,1)
        fn = fns[gids[i]]
        for j in axes(β,2)
            isapprox(β[i,j], 0; atol) && continue
            push!(fn, Shape{T}(cell, lids[j], β[i,j]))
        end
    end
end


function gwpcurl(mesh, edges=nothing; order, interpolator=false)

    T = coordtype(mesh)
    atol = sqrt(eps(T))

    if edges == nothing
        edges = setminus(skeleton(mesh, 1), boundary(mesh))
    end

    C12 = connectivity(mesh, edges, abs)
    R12 = rowvals(C12)
    V12 = nonzeros(C12)

    ne = order+1
    nf = order*(order+1)

    Ne = ne * length(edges)
    Nf = nf * length(mesh)
    Nt = Ne + Nf

    localspace = GWPCurlRefSpace{T,order}()
    refdom = domain(chart(mesh, first(mesh)))
    localdim = numfunctions(localspace, refdom)

    fns = [Shape{T}[] for n in 1:Nt]
    pos = fill(point(0,0,0), Nt)
    for cell in mesh
        cell_ch = chart(mesh, cell)
        for k in nzrange(C12, cell)
            e = R12[k]
            i = V12[k]
            gids = (e-1)*ne .+ (1:ne)
            lids = localindices(localspace, refdom, Val{1}, i)
            edge_ch = chart(edges, e)
            vals = globaldofs(edge_ch, cell_ch, localspace, GWPGlobalEdgeDoFs(order))
            α = vals[lids,:]
            _addshapes!(fns, cell, gids, lids, inv(α'))
        end

        order < 1 && continue
        face_ch = chart(mesh, cell)
        gids = Ne + (cell-1)*nf .+ (1:nf)
        lids = localindices(localspace, refdom, Val{2}, 1)
        vals = globaldofs(face_ch, face_ch, localspace, GWPGlobalFaceDoFs(order))
        α = vals[lids,:]
        _addshapes!(fns, cell, gids, lids, inv(α'))
    end

    for e in edges pos[(e-1)*ne .+ (1:ne)] .= Ref(cartesian(center(chart(edges, e)))) end
    for f in mesh pos[Ne+(f-1)*nf .+ (1: nf)] .= Ref(cartesian(center(chart(mesh, f)))) end

    if interpolator
        return SpaceWithDoFs(
            GWPCurlSpace(mesh, fns, pos, order),
            GWPInterpolator(order, mesh, edges, mesh)
        )
    end
    return GWPCurlSpace(mesh, fns, pos, order)
end

@testitem "GWPcurl global: numfunctions" begin
    using CompScienceMeshes

    h = 0.5
    mesh = meshrectangle(1.0, 1.0, 0.5)
    edges = setminus(skeleton(mesh,1), boundary(mesh))

    order = 2
    gwp = BEAST.gwpcurl(mesh, edges; order=order)

    ne = order+1
    nf = order * (order+1)
    Nt = length(edges)*ne + length(mesh)*nf
    @test numfunctions(gwp) == Nt
end


struct GWPInterpolator
    order
    mesh
    edges
    faces
end


function (interp::GWPInterpolator)(f)

    T = scalartype(f)
    order = interp.order
    mesh = interp.mesh
    edges = interp.edges
    faces = interp.faces

    ne = order+1
    nf = order*(order+1)

    Ne = ne * length(edges)
    Nf = nf * length(faces)
    Nt = Ne + Nf

    vals = zeros(T, Nt)

    C12 = connectivity(edges, mesh, abs)
    R12 = rowvals(C12)
    # V12 = nonzeros(C12)
    for e in edges
        edge_ch = chart(edges, e)
        for k in nzrange(C12, e)
            # l = V12[k]
            c = R12[k]
            cell_ch = chart(mesh, c)
            h = trace(edge_ch, cell_ch, f)
            g = u -> (h(u),)

            gids = (e-1)*ne .+ (1:ne)
            Q = zeros(T, 1, ne)
            globaldofs!(Q, edge_ch, g, GWPGlobalEdgeDoFs(order))
            vals[gids] .= Q[1,:]
            break
        end
    end

    C22 = connectivity(faces, mesh, abs)
    R22 = rowvals(C22)
    for p in faces
        face_ch = chart(faces, p)
        for k in nzrange(C22, p)
            c = R22[k]
            cell_ch = chart(mesh, c)
            h = trace(face_ch, cell_ch, f)
            g = u -> (h(u),)
            gids = Ne + (p-1)*nf .+ (1:nf)
            Q = zeros(T, 1, nf)
            globaldofs!(Q, face_ch, g, GWPGlobalFaceDoFs(order))
            vals[gids] .= Q[1,:]
            break
        end
    end

    return vals
end

# Add to the sparse triplets the result of applying the dofs `Q` (one column per
# dof) to the local shapes on active cell `a`, storing it in rows
# `offset .+ axes(Q,2)`.
function _addinterpentries!(rows, cols, vals, ad, a, Q, offset)
    for j in axes(Q,2)
        for i in axes(Q,1)
            for (m,α) in ad[a,i]
                iszero(α) && continue
                push!(rows, offset + j)
                push!(cols, m)
                push!(vals, α * Q[i,j])
            end
        end
    end
end

# Return the active index of the first cell adjacent to column `k` of `C` that
# supports at least one function, or 0 if there is no such cell.
function _firstactivecell(C, k, g2a)
    R = rowvals(C)
    for m in nzrange(C, k)
        a = get(g2a, R[m], 0)
        a != 0 && return a
    end
    return 0
end

"""
    (interp::GWPInterpolator)(space::Space)

Interpolate each of the generating functions in `space` on a GWP space. It
is assumed that the GWP space geometry is defined on the same vertex buffer
as `space`, and that `space` is tangentially continuous, so that the edge dofs
can be evaluated from either side of an edge.

Returns a sparse matrix `A` of size `Nt × numfunctions(space)`, where column `m`
contains the GWP coefficients of the interpolant of the `m`-th function.
"""
function (interp::GWPInterpolator)(space::Space)

    T = scalartype(space)
    edges = interp.edges
    faces = interp.faces

    order = interp.order
    ne = order+1
    nf = order*(order+1)

    Ne = ne * length(edges)
    Nf = nf * length(faces)
    Nt = Ne + Nf

    rows = Int[]
    cols = Int[]
    vals = T[]

    charts1, ad, a2g = assemblydata(space)
    g2a = Dict(p => i for (i,p) in enumerate(a2g))
    ϕ = refspace(space)

    cells1 = geometry(space)
    D12 = connectivity(edges, cells1, abs)
    D22 = connectivity(faces, cells1, abs)

    for e in edges
        # functions of `space` vanish on cells where none is supported
        a = _firstactivecell(D12, e, g2a)
        a == 0 && continue

        edge_ch = chart(edges, e)
        Q = globaldofs(edge_ch, charts1[a], ϕ, GWPGlobalEdgeDoFs(order))
        _addinterpentries!(rows, cols, vals, ad, a, Q, ne*(e-1))
    end

    if order > 0
        for p in faces
            a = _firstactivecell(D22, p, g2a)
            a == 0 && continue

            face_ch = chart(faces, p)
            Q = globaldofs(face_ch, charts1[a], ϕ, GWPGlobalFaceDoFs(order))
            _addinterpentries!(rows, cols, vals, ad, a, Q, Ne + nf*(p-1))
        end
    end

    return sparse(rows, cols, vals, Nt, numfunctions(space))
end
