abstract type HHHFunctional{T} <: Functional end

mutable struct FunctionExcitation{T} <: HHHFunctional{T}
    f::Function
    coeff::T
end
mutable struct NdotExcitation{T} <: HHHFunctional{T}
    e::HHHFunctional{T}
end
mutable struct NcrossExcitation{T} <: HHHFunctional{T}
    e::HHHFunctional{T}
end

FunctionExcitation(f) = FunctionExcitation(f,f([0.0,0.0,0.0]))

function (e::FunctionExcitation)(p)
    x = cartesian(p)
    e.f(x)
end
function (e::NdotExcitation)(p)
    n = normal(p)
    dot(n,e.e(p))
end
function (e::NcrossExcitation)(p)
    n = normal(p)
    cross(n,e.e(p))
end



integrand(::HHHFunctional,test_vals, field_val) = dot(test_vals.value,field_val)


scalartype(e::HHHFunctional{T}) where {T<:Number} = T
scalartype(e::HHHFunctional{Vector{T}}) where {T<:Number} = T