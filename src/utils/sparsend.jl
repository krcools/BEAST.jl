module SparseND

struct Banded3D{T} <: AbstractArray{T,3}
    k0::Array{Int,2}
    data::Array{T,3}
    maxk0::Int
    maxk1::Int
    function Banded3D{T}(k0,data,maxk1) where T
        @assert size(k0) == size(data)[2:3]
        return new(k0,data,maximum(k0),maxk1)
    end
end

Banded3D(k0,data::Array{T},maxk1) where {T} = Banded3D{T}(k0,data,maxk1)

bandwidth(A::Banded3D) = size(A.data,1)

import Base: size, getindex, setindex!

# size(A::Banded3D) = size(A.data)
size(A::Banded3D) = tuple(size(A.k0)..., A.maxk1)

function getindex(A::Banded3D, m::Int, n::Int, k::Int)
    k0 = A.k0[m,n]
    l = k - k0 + 1
    l < 1 && return zero(eltype(A))
    l > bandwidth(A) && return zero(eltype(A))
    A.data[l,m,n]
end

function setindex!(A::Banded3D, v, m, n, k)
    @assert A.k0[m,n] <= k <= A.k0[m,n] + bandwidth(A) - 1
    A.data[k-A.k0[m,n]+1,m,n] = v
end

# """
#     convolve{T}(A::Banded3D, x::Array{T,2})
#
# Compute the *space-time* convolution of A and x.
#
# Returns an array `y` of size `(size(A,1),size(x,2))` such that
# `y[m,i] = sum([A[]])``
# """
# function convolve(A::Banded3D, x::Array{T,2}) where T
#     V = promote_type(eltype(A),eltype(x))
#     y = zeros(V,size(A,1),size(x,2))
#     bw = bandwidth(A)
#     for i in 1:size(A,1)
#         for l in 1 : size(y,2)
#             for j in 1:size(A,2)
#                 k0 = A.k0[i,j]
#                 for k in k0 : k0 + bw - 1
#                     l-k+1 < 1 && continue
#                     y[i,l] += A.data[k,i,j] * x[j,l-k+1]
#                 end
#             end
#         end
#     end
#     return y
# end

end # module
