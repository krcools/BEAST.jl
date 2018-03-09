"""
Array type of rank 3 (taking three scalar indices) that is used for the storage
of discrete convolution operators that happen to be the tensor product of a spatial
part and a temporal part.

Using a dedicated storage type for these operators saves both on memory and computation
time for convolutions with space-time data.
"""
mutable struct SpaceTimeDyad{A,B,T} <: AbstractArray{3,T}
    spatial_factor::A
    temporal_factor::B
end

# Implementation of a minimalistic array interface
# Note that setindex! is not supported. There is no way to enforce maintaining
# the dyadic structure during element modification.
import Base: size, getindex, linearindexing
size(X::SpaceTimeDyad) = (size(X.spatial_factor)..., length(X.temporal_factor))
linearindexing(::Type{T}) where {T<:SpaceTimeDyad} = Base.LinearSlow()
function getindex(X::SpaceTimeDyad, m::Int, n::Int, k::Int)
    T = promote_type(eltype(X.spatial_factor), eltype(X.temporal_factor))
    k > length(X.spacetial_factor) && return zero(T)
    return X.spatial_factor[m,n] * X.temporal_factor[k]
end


function convolve(X::SpaceTimeDyad, x::Array{T,2})
    
end
