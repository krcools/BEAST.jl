struct FunctionExcitation{T} <: Functional
    f
end
export FunctionExcitation
kernelvals(op::ComposedOperator,a) = nothing



function (func::FunctionExcitation)(x)
    return func.f(x)
end
function (func::FunctionExcitation)(x::Union{CompScienceMeshes.MeshPointNM,BEAST.TraceMeshPointNM})
    return func.f(cartesian(x))
end
scalartype(ff::FunctionExcitation{T}) where {T} = T
cross(::NormalVector, p::FunctionExcitation) = CrossTraceMW(p)
integrand(::FunctionExcitation,tval,fval) = tval[1]*fval
export FunctionExcitation


function farfieldlocal!(zlocal,op::ComposedOperatorIntegral,refspace,y,el,qr)

    for q in qr
        x = q.point
        F = q.value
        dx = q.weight

        krn = kernelvals(op, y, x)
        int = integrand(op,krn,y,F,x) * dx
        for (i,j) in enumerate(int)
            zlocal[i] += j
        end

    end

end