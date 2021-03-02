module SparseND

import ...BEAST

struct MatrixOfConvolutions{T} <: AbstractArray{Vector{T},2}
    # banded::AbstractArray{T,3}
    convop::BEAST.ConvOp{T}
end

function Base.eltype(x::MatrixOfConvolutions{T}) where {T}
    Vector{T}
end
Base.size(x::MatrixOfConvolutions) = size(x.convop)[1:2]
function Base.getindex(x::MatrixOfConvolutions, m, n)
    # return x.banded[m,n,:]
    return x.convop[m,n,:]
end


struct SpaceTimeData{T} <: AbstractArray{Vector{T},1}
    data::Array{T,2}
end

Base.eltype(x::SpaceTimeData{T}) where {T} = Vector{T}
Base.size(x::SpaceTimeData) = (size(x.data)[1],)
Base.getindex(x::SpaceTimeData, i::Int) = x.data[i,:]

end # module
