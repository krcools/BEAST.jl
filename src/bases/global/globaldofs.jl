function trace(edge_ch, face_ch, fields)
    T = coordtype(edge_ch)
    atol = sqrt(eps(T))

    face_dom = domain(face_ch)

    edge_verts = vertices(edge_ch)
    face_verts = vertices(face_ch)

    perm = [something(findfirst(w -> isapprox(v, w; atol), face_verts), 0) for v in edge_verts]
    # @show perm
    injection = simplex(vertices(face_dom)[perm])

    function f(s)
        u = neighborhood(injection, s)
        p = neighborhood(face_ch, cartesian(u))
        fields(p)
    end

    return f
end



struct _LagrangeGlobalEdgeDoFs
    order::Int
end
function numfunctions(dof::_LagrangeGlobalEdgeDoFs) dof.order-1 end
function (dof::_LagrangeGlobalEdgeDoFs)(s)
    T = typeof(s)
    nodes = range(zero(T), one(T), length=order+1)
    [_lagpoly(nodes, i, s) for i in 2:order]
end

struct _LagrangeGlobalFaceDoFs
    order::Int
end
# function numfunctions(dof::_LagrangeGlobalFaceDoFs) div((order-1)*(order-2),2) end
# function (dof::_LagrangeGlobalFaceDoFs)(s)
#     T = eltype(s)
#     nodes = range(zero(T), one(T), length=order+1)
#     r = zeros(T, numfunctions(dof))
#     s1, s2 = s
#     s3 = 1 - s1 - s2
#     idx = 1
#     degree = dof.degree
#     for i in 0:degree
#         prodi = _lagpoly(nodes, i+1, s1, i)
#         for j in 0:degree
#             k = degree - i - j
#             k < 0 && continue
#             prodj = _lagpoly(nodes, j+1, s2, j)
#             prodk = _lagpoly(nodes, k+1, s3, k)
#             r[idx] = prodi * prodj * prodk
#     end end
#     return r
# end

function globaldofs(edge_ch, face_ch, localspace, dof::_LagrangeGlobalEdgeDoFs)

    T = coordtype(edge_ch)
    f = trace(edge_ch, face_ch, localspace)
    # r = zero(numfunctions(localspace), numfunctions(dof))
    # for (s,w) in CompScienceMeshes.quadpoints(edge_dom, 2*order)
    #     u = neighborhood(injection, s)
    #     p = neighborhood(face_ch, cartesian(u))
    #     r .+= w * [x.value for x in localspace(p)] * dof(s)'
    # end
    
    ds = one(T) / dof.order
    return stack(range(ds, step=ds, length=dof.order-1)) do s
        [x.value for x in f(s)]
    end
end


function globaldofs(edge_ch, face_ch, localspace, dof::_LagrangeGlobalFaceDoFs)

    T = coordtype(edge_ch)
    f = trace(edge_ch, face_ch, localspace)
    # r = zero(numfunctions(localspace), numfunctions(dof))
    # for (s,w) in CompScienceMeshes.quadpoints(edge_dom, 2*order)
    #     u = neighborhood(injection, s)
    #     p = neighborhood(face_ch, cartesian(u))
    #     r .+= w * [x.value for x in localspace(p)] * dof(s)'
    # end
    
    d = dof.order
    S = ((i,j,d-i-j) for i in 0:d for j in 0:d if (i+j < d && i > 0 && j > 0))
    return stack(S) do s
        s = (s[1]/d,s[2]/d)
        [x.value for x in f(s)]
    end
end

struct _LagrangeGlobalNodesDoFs
    order::Int
end

function globaldofs(edge_ch, face_ch, localspace, dof::_LagrangeGlobalNodesDoFs)
    f = trace(edge_ch, face_ch, localspace)
    return stack([()]) do s
        [x.value for x in f(s)]
    end
end


@testitem "globaldofs: interpolatory on edge" begin
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

    localspace = BEAST.LagrangeRefSpace{T,3,3,10}()
    L = BEAST._LagrangeGlobalEdgeDoFs(3)
    dofs = BEAST.globaldofs(edge_ch, face_ch, localspace, L)

    @test size(dofs) == (10,2)
    A = [
        0 0
        0 1
        1 0
        0 0
        0 0
        0 0
        0 0
        0 0
        0 0
        0 0]
    @test dofs ≈ A
end