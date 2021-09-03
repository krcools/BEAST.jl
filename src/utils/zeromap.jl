import LinearMaps

struct ZeroMap{T} <: LinearMap{T}
    num_rows::Int
    num_cols::Int
end

Base.size(A::ZeroMap) = (A.num_rows, A.num_cols,)

function LinearAlgebra.mul!(y::AbstractVector, L::ZeroMap, x::AbstractVector, α::Number, β::Number)
    y .= β * y
    # LinearAlgebra.mul!(yI, AIJ, xJ, α, β)
end