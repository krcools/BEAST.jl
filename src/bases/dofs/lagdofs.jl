struct _LagrangeGlobalEdgeDoFs
    order::Int
end
function numfunctions(dof::_LagrangeGlobalEdgeDoFs) dof.order-1 end
function (dof::_LagrangeGlobalEdgeDoFs)(s)
    T = typeof(s)
    order = dof.order
    nodes = range(zero(T), one(T), length=order+1)
    [_lagpoly(nodes, i, s) for i in 2:order]
end

struct _LagrangeGlobalFaceDoFs
    order::Int
end
function numfunctions(dof::_LagrangeGlobalFaceDoFs) div((dof.order-1)*(dof.order-2), 2) end

struct _LagrangeGlobalNodesDoFs
    order::Int
end
function numfunctions(dof::_LagrangeGlobalNodesDoFs) 1 end


"""
    globaldofs(global_ch, face_ch, localspace, dof)

Evaluate the global Lagrange degrees of freedom `dof` defined on `global_ch` (a
vertex, edge, or face) on the local basis functions in `localspace` defined on
`face_ch`. The returned matrix has one row per local basis function and one
column per degree of freedom.
"""
function globaldofs(global_ch, face_ch, localspace, dof::Union{
        _LagrangeGlobalNodesDoFs, _LagrangeGlobalEdgeDoFs, _LagrangeGlobalFaceDoFs})

    T = coordtype(global_ch)
    f = trace(global_ch, face_ch, localspace)
    g = u -> collect(fi.value for fi in f(u))
    Q = zeros(T, numfunctions(localspace, domain(face_ch)), numfunctions(dof))
    return globaldofs!(Q, global_ch, g, dof)
end

# The DoFs are point evaluations. `pulled_back_fn(u)` returns the values at
# the point with parameters `u` on `global_ch` of the functions corresponding
# to the rows of `Q`.
function globaldofs!(Q, global_ch, pulled_back_fn, dof::_LagrangeGlobalNodesDoFs)
    vals = pulled_back_fn(())
    for j in axes(Q,1)
        Q[j,1] = vals[j]
    end
    return Q
end

function globaldofs!(Q, global_ch, pulled_back_fn, dof::_LagrangeGlobalEdgeDoFs)

    d = dof.order
    T = coordtype(global_ch)
    ds = one(T) / d

    for i in 1:d-1
        vals = pulled_back_fn(i*ds)
        for j in axes(Q,1)
            Q[j,i] = vals[j]
        end
    end
    return Q
end

function globaldofs!(Q, global_ch, pulled_back_fn, dof::_LagrangeGlobalFaceDoFs)

    d = dof.order
    T = coordtype(global_ch)

    S = ((i,j) for i in 0:d for j in 0:d if (i+j < d && i > 0 && j > 0))
    for (idx, (i,j)) in enumerate(S)
        vals = pulled_back_fn((T(i)/d, T(j)/d))
        for k in axes(Q,1)
            Q[k,idx] = vals[k]
        end
    end
    return Q
end


"""
    LagrangeInterpolator(order, mesh, nodes, edges, faces)

Interpolator onto the continuous Lagrange space of degree `order` built by
`lagrangec0(mesh, nodes, edges; order)`. Calling it on a scalar field `f`
returns the coefficients of the interpolant, i.e. the values of the vertex,
edge, and face based degrees of freedom applied to `f`, in the same order as
the global basis functions. Normally `faces` is `mesh`.
"""
struct LagrangeInterpolator
    order
    mesh
    nodes
    edges
    faces
end


function (interp::LagrangeInterpolator)(f)

    T = scalartype(f)
    order = interp.order
    mesh = interp.mesh
    nodes = interp.nodes
    edges = interp.edges
    faces = interp.faces

    ne = order-1
    nf = div((order-1)*(order-2), 2)

    nV = length(nodes)
    nE = ne * length(edges)
    nF = nf * length(faces)

    vals = zeros(T, nV + nE + nF)

    C02 = connectivity(nodes, mesh, abs)
    R02 = rowvals(C02)
    for v in nodes
        node_ch = chart(nodes, v)
        for k in nzrange(C02, v)
            c = R02[k]
            cell_ch = chart(mesh, c)
            h = trace(node_ch, cell_ch, f)
            g = u -> (h(u),)
            Q = zeros(T, 1, 1)
            globaldofs!(Q, node_ch, g, _LagrangeGlobalNodesDoFs(order))
            vals[v] = Q[1,1]
            break
        end
    end

    order < 2 && return vals
    C12 = connectivity(edges, mesh, abs)
    R12 = rowvals(C12)
    for e in edges
        edge_ch = chart(edges, e)
        for k in nzrange(C12, e)
            c = R12[k]
            cell_ch = chart(mesh, c)
            h = trace(edge_ch, cell_ch, f)
            g = u -> (h(u),)
            gids = nV + (e-1)*ne .+ (1:ne)
            Q = zeros(T, 1, ne)
            globaldofs!(Q, edge_ch, g, _LagrangeGlobalEdgeDoFs(order))
            vals[gids] .= Q[1,:]
            break
        end
    end

    order < 3 && return vals
    C22 = connectivity(faces, mesh, abs)
    R22 = rowvals(C22)
    for p in faces
        face_ch = chart(faces, p)
        for k in nzrange(C22, p)
            c = R22[k]
            cell_ch = chart(mesh, c)
            h = trace(face_ch, cell_ch, f)
            g = u -> (h(u),)
            gids = nV + nE + (p-1)*nf .+ (1:nf)
            Q = zeros(T, 1, nf)
            globaldofs!(Q, face_ch, g, _LagrangeGlobalFaceDoFs(order))
            vals[gids] .= Q[1,:]
            break
        end
    end

    return vals
end


