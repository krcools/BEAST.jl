import LinearMaps

struct ZeroMap{T,U,V} <: LinearMap{T}
    range::U
    domain::V
end

ZeroMap{T}(range::U, domain::V) where {T,U,V} = ZeroMap{T,U,V}(range, domain)
LinearMaps.MulStyle(A::ZeroMap) = LinearMaps.FiveArg()

Base.size(A::ZeroMap) = (length(A.range), length(A.domain),)
Base.axes(A::ZeroMap) = (A.range, A.domain)

function LinearMaps._unsafe_mul!(y::AbstractVector, L::ZeroMap, x::AbstractVector, α::Number, β::Number)
    y .*= β
end

# function LinearAlgebra.mul!(Y::PseudoBlockMatrix, c::Number, X::ZeroMap, a::Number, b::Number)
#     @assert b == 1
#     return Y
# end

function LinearMaps._unsafe_mul!(Y::AbstractMatrix, X::ZeroMap, c::Number, a::Number, b::Number)
    rmul!(Y, b)
    return Y
end

function LinearMaps._unsafe_mul!(Y::AbstractMatrix, X::ZeroMap, c::Number)
    fill!(Y, false)
    return Y
end

# function LinearAlgebra.mul!(y::AbstractVector,
#     Lt::LinearMaps.TransposeMap{<:Any,<:ZeroMap},
#     x::AbstractVector, α::Number, β::Number)

#     y .*= β
# end


function LinearMaps._unsafe_mul!(y::AbstractVector, L::ZeroMap, x::AbstractVector)
    y .= 0
end

# function LinearAlgebra.mul!(y::AbstractVector,
#     Lt::LinearMaps.TransposeMap{<:Any,<:ZeroMap}, x::AbstractVector)

#     y .= 0
# end

LinearAlgebra.adjoint(A::ZeroMap{T}) where {T} = ZeroMap{T}(A.domain, A.range)
LinearAlgebra.transpose(A::ZeroMap{T}) where {T} = ZeroMap{T}(A.domain, A.range)

# function Base.:(*)(A::ZeroMap, x::AbstractVector)
#     T = eltype(A)
#     y = similar(A.range, T)
#     fill!(y, zero(T))
#     LinearAlgebra.mul!(y,A,x)
# end

# Base.Matrix{T}(A::ZeroMap) where {T} = zeros(T, length(A.range), length(A.domain))
