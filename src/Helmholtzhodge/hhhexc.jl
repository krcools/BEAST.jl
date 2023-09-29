mutable struct FunctionExcitation{T} <: Functional
    f::Function
    coeff::T
end
FunctionExcitation(f) = FunctionExcitation(f,f([0.0,0.0,0.0]))

function (e::FunctionExcitation)(p)
    x = cartesian(p)
    e.f(x)
end
integrand(::FunctionExcitation,test_vals, field_val) = dot(test_vals.value,field_val)
scalartype(e::FunctionExcitation{T}) where {T<:Number} = T
scalartype(e::FunctionExcitation{Vector{T}}) where {T<:Number} = T