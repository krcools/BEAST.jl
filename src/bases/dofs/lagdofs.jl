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

function globaldofs(edge_ch, face_ch, localspace, dof::_LagrangeGlobalEdgeDoFs)

    T = coordtype(edge_ch)
    f = trace(edge_ch, face_ch, localspace)
    
    ds = one(T) / dof.order
    return stack(range(ds, step=ds, length=dof.order-1)) do s
        [x.value for x in f(s)]
    end
end


function globaldofs(edge_ch, face_ch, localspace, dof::_LagrangeGlobalFaceDoFs)

    T = coordtype(edge_ch)
    f = trace(edge_ch, face_ch, localspace)
    
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