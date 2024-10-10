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

function globaldofs(edge_ch, face_ch, localspace, dof::GWPGlobalEdgeDoFs)
    f = trace(edge_ch, face_ch, localspace)
    T = coordtype(edge_ch)
    ds = one(T) / (dof.order+2)
    return stack(range(ds, step=ds, length=dof.order+1)) do s
        p = neighborhood(edge_ch, s)
        t = tangents(p, 1)
        [dot(t,x.value) for x in f(s)]
    end
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
    display(round.(A, digits=3))
end

function globaldofs(edge_ch, face_ch, localspace, dof::GWPGlobalFaceDoFs)

    d = dof.order
    T = coordtype(edge_ch)
    ds = one(T) / (d+2)

    f = trace(edge_ch, face_ch, localspace)
    Q = zeros(T, numfunctions(localspace, domain(face_ch)), d*(d+1))
    S = ((i,j,d+2-i-j) for i in 1:d+1 for j in 1:d+1 if d+2-i-j > 0)
    # return stack(S) do (i,j,k)
    idx = 1
    for (i,j,k) in S
        u = (i*ds,j*ds)
        p = neighborhood(edge_ch, u)
        t1 = tangents(p,1)
        t2 = tangents(p,2)
        vals = f(u)
        q1 = [dot(t1,x.value) for x in vals]
        q2 = [dot(t2,x.value) for x in vals]
        Q[:,idx] = q1; idx += 1
        Q[:,idx] = q2; idx += 1
    end

    return Q
end


@testitem "GWP face global dofs" begin
    using CompScienceMeshes
    T = Float64
    face_ch = simplex(
        point(3,0,0),
        point(2,0,-1),
        point(0,0,1),
    )
    # edge_ch = simplex(
    #     point(0,0,1),
    #     point(2,0,-1),
    #     point(3,0,0)
    # )
    edge_ch = simplex(
        point(3,0,0),
        point(2,0,-1),
        point(0,0,1),
    )

    degree = 3
    localspace = BEAST.GWPCurlRefSpace{T,degree}()
    dofs = BEAST.GWPGlobalFaceDoFs(degree)
    A = BEAST.globaldofs(edge_ch, face_ch, localspace, dofs)
    display(round.(A, digits=3))
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


function gwpcurl(mesh, edges=nothing; order)

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

        order < 2 && continue
        face_ch = chart(mesh, cell)
        gids = Ne + (cell-1)*nf .+ (1:nf)
        lids = localindices(localspace, refdom, Val{2}, 1)
        vals = globaldofs(face_ch, face_ch, localspace, GWPGlobalFaceDoFs(order))
        α = vals[lids,:]
        _addshapes!(fns, cell, gids, lids, inv(α'))
    end

    for e in edges pos[(e-1)*ne .+ (1:ne)] .= Ref(cartesian(center(chart(edges, e)))) end
    for f in mesh pos[Ne+(f-1)*nf .+ (1: nf)] .= Ref(cartesian(center(chart(mesh, f)))) end

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