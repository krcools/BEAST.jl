import LinearMaps

struct Rank1Map{T,U,V} <: LinearMap{T}
    u::U
    v::V
end

Rank1Map{T}(u::U, v::V) where {T,U,V} = Rank1Map{T,U,V}(u,v)
LinearMaps.MulStyle(A::Rank1Map) = LinearMaps.FiveArg()

Base.size(A::Rank1Map) = (length(A.u), length(A.v),)
Base.axes(A::Rank1Map) = (axes(A.u)..., axes(A.v)...)

function LinearMaps._unsafe_mul!(y::AbstractVector, L::Rank1Map, x::AbstractVector, α::Number=true, β::Number=false)
    y .*= β
    y .+= α .* L.u .* dot(L.v, x)
end

function LinearMaps._unsafe_mul!(Y::AbstractMatrix, L::Rank1Map, c::Number, a::Number=true, b::Number=false)
    rmul!(Y, b)
    Y .+= (L.u * L.v') .* c .* α
    return Y
end
